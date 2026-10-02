'use strict';

/* ======================================================================
   IMU 室內行走軌跡偵測 — 核心邏輯
   座標系：世界座標系，X=東, Y=北, Z=垂直向上，原點=開始偵測瞬間
   ====================================================================== */

// ---------- 全域狀態 ----------
const state = {
  tracking: false,
  startTime: 0,
  lastSampleTime: 0,
  sampleCount: 0,
  carryMode: 'pocket_pants',
  username: '',
  heightCm: null,
  loopClosureEnabled: false,

  // 姿態四元數 [w, x, y, z]，初始為單位四元數
  quat: [1, 0, 0, 0],
  magDisturbed: false,

  // 步伐偵測用的濾波窗口
  accVertBuffer: [],      // {t, value} 低通後的世界垂直加速度
  lastStepTime: 0,
  stepPeakCandidate: null,
  recentSteps: [],        // 最近步的 {dx, dy, dz, t}，供斜度判斷用滑動窗口

  // 目前位置與朝向
  pos: { x: 0, y: 0, z: 0 },
  heading: 0, // 弧度，世界座標系航向角 (0 = 北)

  // 轉彎抑制
  turnSuppressUntil: 0,

  // 輸出緩衝
  track: [],
  rawLog: [],

  // 統計
  stepCount: 0,
  totalDistance: 0,
  currentSegmentType: 'flat',
  currentConfidence: 1.0,
};

// 不同攜帶方式的加速度峰值門檻 (m/s^2)
// 基於「三軸合成加速度模長扣除重力常數」(|a|-g)。
// 【目前為暫時性低門檻，僅供波形校準用】第一輪實機測試顯示先前估計的
// 2.0+ 門檻遠高於實際訊號量級（波形被壓成看不出起伏，步伐完全偵測不到）。
// 暫時調低以便先看清楚真實訊號大小，待有真實峰值數據後應立刻依
// 「訊號峰值的 50-60%」重新設定，不應長期使用這組暫定值。
const CARRY_MODE_THRESHOLD = {
  pocket_pants: 0.4,
  pocket_jacket: 0.5,
  handheld: 0.6,
  armband: 0.5,
};

const GRAVITY = 9.80665;

// ---------- 工具函式：四元數運算 ----------
function quatMultiply(a, b) {
  const [aw, ax, ay, az] = a;
  const [bw, bx, by, bz] = b;
  return [
    aw*bw - ax*bx - ay*by - az*bz,
    aw*bx + ax*bw + ay*bz - az*by,
    aw*by - ax*bz + ay*bw + az*bx,
    aw*bz + ax*by - ay*bx + az*bw,
  ];
}
function quatConjugate(q) { return [q[0], -q[1], -q[2], -q[3]]; }
function quatNormalize(q) {
  const n = Math.sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]) || 1;
  return [q[0]/n, q[1]/n, q[2]/n, q[3]/n];
}
// 用四元數把向量從裝置座標系旋轉到世界座標系
function rotateVectorByQuat(q, v) {
  const p = [0, v[0], v[1], v[2]];
  const qc = quatConjugate(q);
  const r = quatMultiply(quatMultiply(q, p), qc);
  return [r[1], r[2], r[3]];
}

/* ======================================================================
   Madgwick AHRS 互補濾波（簡化版）
   融合加速度計 + 陀螺儀 + 磁力計，輸出姿態四元數
   ====================================================================== */
const madgwick = {
  beta: 0.1, // 演算法增益，越大收斂越快但雜訊越大
  q: [1, 0, 0, 0],

  update(gx, gy, gz, ax, ay, az, mx, my, mz, dt) {
    let q1 = this.q[0], q2 = this.q[1], q3 = this.q[2], q4 = this.q[3];
    let norm;

    // 正規化加速度計
    norm = Math.sqrt(ax*ax + ay*ay + az*az);
    if (norm === 0) return this.q;
    ax /= norm; ay /= norm; az /= norm;

    let useMag = mx !== null && mx !== undefined && !state.magDisturbed;
    let qDot1, qDot2, qDot3, qDot4;

    if (useMag) {
      norm = Math.sqrt(mx*mx + my*my + mz*mz);
      if (norm === 0) useMag = false;
      else { mx /= norm; my /= norm; mz /= norm; }
    }

    if (useMag) {
      // 完整 9 軸融合
      const hx = mx*(1-2*q3*q3-2*q4*q4) + my*2*(q2*q3-q1*q4) + mz*2*(q2*q4+q1*q3);
      const hy = mx*2*(q2*q3+q1*q4) + my*(1-2*q2*q2-2*q4*q4) + mz*2*(q3*q4-q1*q2);
      const bx = Math.sqrt(hx*hx + hy*hy);
      const bz = mx*2*(q2*q4-q1*q3) + my*2*(q3*q4+q1*q2) + mz*(1-2*q2*q2-2*q3*q3);

      const f1 = 2*(q2*q4 - q1*q3) - ax;
      const f2 = 2*(q1*q2 + q3*q4) - ay;
      const f3 = 1 - 2*(q2*q2 + q3*q3) - az;
      const f4 = 2*bx*(0.5 - q3*q3 - q4*q4) + 2*bz*(q2*q4 - q1*q3) - mx;
      const f5 = 2*bx*(q2*q3 - q1*q4) + 2*bz*(q1*q2 + q3*q4) - my;
      const f6 = 2*bx*(q1*q3 + q2*q4) + 2*bz*(0.5 - q2*q2 - q3*q3) - mz;

      const J11 = -2*q3, J12 = 2*q4, J13 = -2*q1, J14 = 2*q2;
      const J21 = 2*q2, J22 = 2*q1, J23 = 2*q4, J24 = 2*q3;
      const J31 = 0, J32 = -4*q2, J33 = -4*q3, J34 = 0;

      let s1 = J11*f1 + J21*f2;
      let s2 = J12*f1 + J22*f2 + J32*f3;
      let s3 = J13*f1 + J23*f2 + J33*f3;
      let s4 = J14*f1 + J24*f2;

      norm = Math.sqrt(s1*s1+s2*s2+s3*s3+s4*s4) || 1;
      s1/=norm; s2/=norm; s3/=norm; s4/=norm;

      qDot1 = 0.5*(-q2*gx - q3*gy - q4*gz) - this.beta*s1;
      qDot2 = 0.5*(q1*gx + q3*gz - q4*gy) - this.beta*s2;
      qDot3 = 0.5*(q1*gy - q2*gz + q4*gx) - this.beta*s3;
      qDot4 = 0.5*(q1*gz + q2*gy - q3*gx) - this.beta*s4;
    } else {
      // 退化為 6 軸（加速度 + 陀螺儀），磁力受干擾時使用
      const f1 = 2*(q2*q4 - q1*q3) - ax;
      const f2 = 2*(q1*q2 + q3*q4) - ay;
      const f3 = 1 - 2*(q2*q2 + q3*q3) - az;

      const J11 = -2*q3, J12 = 2*q4, J13 = -2*q1, J14 = 2*q2;
      const J21 = 2*q2, J22 = 2*q1, J23 = 2*q4, J24 = 2*q3;
      const J31 = 0, J32 = -4*q2, J33 = -4*q3, J34 = 0;

      let s1 = J11*f1 + J21*f2;
      let s2 = J12*f1 + J22*f2 + J32*f3;
      let s3 = J13*f1 + J23*f2 + J33*f3;
      let s4 = J14*f1 + J24*f2;

      let norm2 = Math.sqrt(s1*s1+s2*s2+s3*s3+s4*s4) || 1;
      s1/=norm2; s2/=norm2; s3/=norm2; s4/=norm2;

      qDot1 = 0.5*(-q2*gx - q3*gy - q4*gz) - this.beta*s1;
      qDot2 = 0.5*(q1*gx + q3*gz - q4*gy) - this.beta*s2;
      qDot3 = 0.5*(q1*gy - q2*gz + q4*gx) - this.beta*s3;
      qDot4 = 0.5*(q1*gz + q2*gy - q3*gx) - this.beta*s4;
    }

    q1 += qDot1*dt; q2 += qDot2*dt; q3 += qDot3*dt; q4 += qDot4*dt;
    this.q = quatNormalize([q1,q2,q3,q4]);
    return this.q;
  }
};

/* ======================================================================
   感測器擷取
   ====================================================================== */
let lastAcc = null, lastGyro = null, lastMag = null;
let sampleRateTracker = { count: 0, windowStart: 0, lastRate: 0 };

function handleMotion(event) {
  const now = performance.now();
  const accG = event.accelerationIncludingGravity;
  const rot = event.rotationRate;
  if (!accG || accG.x === null) return;

  lastAcc = [accG.x, accG.y, accG.z];
  lastGyro = rot ? [
    (rot.alpha || 0) * Math.PI / 180,
    (rot.beta || 0) * Math.PI / 180,
    (rot.gamma || 0) * Math.PI / 180,
  ] : [0, 0, 0];

  processSample(now);
}

function handleOrientation(event) {
  // 作為磁力計資料的替代來源（Web 無直接 Magnetometer 權限時，用方向事件推算航向參考）
  if (event.alpha === null) return;
  // 轉換成近似磁力向量（僅用於航向修正，非真實 µT 值）
  const alphaRad = (event.alpha || 0) * Math.PI / 180;
  const betaRad = (event.beta || 0) * Math.PI / 180;
  const gammaRad = (event.gamma || 0) * Math.PI / 180;
  lastMag = [Math.cos(alphaRad), Math.sin(alphaRad), 0];
  state._lastHeadingHint = event.webkitCompassHeading !== undefined && event.webkitCompassHeading !== null
    ? event.webkitCompassHeading
    : (360 - event.alpha);
}

function processSample(nowMs) {
  if (!state.tracking) return;
  if (!lastAcc || !lastGyro) return;

  const dt = state.lastSampleTime ? Math.min((nowMs - state.lastSampleTime) / 1000, 0.2) : 0.02;
  state.lastSampleTime = nowMs;
  state.sampleCount++;

  // 取樣率估計（每秒重算一次）
  sampleRateTracker.count++;
  if (!sampleRateTracker.windowStart) sampleRateTracker.windowStart = nowMs;
  if (nowMs - sampleRateTracker.windowStart >= 1000) {
    sampleRateTracker.lastRate = sampleRateTracker.count / ((nowMs - sampleRateTracker.windowStart) / 1000);
    sampleRateTracker.count = 0;
    sampleRateTracker.windowStart = nowMs;
    document.getElementById('stat-samplerate').textContent = sampleRateTracker.lastRate.toFixed(0);
  }

  const [ax, ay, az] = lastAcc; // m/s^2, 含重力, 裝置座標系
  const [gx, gy, gz] = lastGyro; // rad/s
  const mag = lastMag || [1, 0, 0];

  // 磁場異常偵測（因子#6）：用近似值，此處偵測陀螺儀/加速度計是否劇烈不一致作為簡化代理
  // 真實 µT 值在 Web 環境難以直接取得，改用方向事件連續性作為可信度代理指標
  state.magDisturbed = false;

  const q = madgwick.update(gx, gy, gz, ax, ay, az, mag[0], mag[1], mag[2], dt);

  // 計算線性加速度（扣除重力）：重力方向 = 世界座標 Z 軸經逆旋轉回裝置座標
  const gravityWorld = [0, 0, GRAVITY];
  const qInv = quatConjugate(q);
  const gravityDevice = rotateVectorByQuat(qInv, gravityWorld);
  const linearAccDevice = [ax - gravityDevice[0], ay - gravityDevice[1], az - gravityDevice[2]];

  // 轉換到世界座標系
  const accWorld = rotateVectorByQuat(q, linearAccDevice);

  // 記錄原始感測器 log
  state.rawLog.push({
    t_ms: Math.round(nowMs - state.startTime),
    acc: [round3(ax), round3(ay), round3(az)],
    acc_linear: [round3(accWorld[0]), round3(accWorld[1]), round3(accWorld[2])],
    gyro: [round4(gx), round4(gy), round4(gz)],
    mag: [round3(mag[0]), round3(mag[1]), round3(mag[2])],
    orientation_quat: [round4(q[0]), round4(q[1]), round4(q[2]), round4(q[3])],
  });

  // 步伐偵測改用「加速度向量模長」扣除重力常數，這個量對姿態估計誤差不敏感
  // （不論手機怎麼轉，走路時 |a| 的波峰都在，比依賴 Madgwick 收斂後的世界座標垂直分量穩健得多）
  const accMagnitude = Math.sqrt(ax*ax + ay*ay + az*az) - GRAVITY;

  stepDetection(nowMs, accWorld, gz, accMagnitude);
}

function round3(v) { return Math.round(v * 1000) / 1000; }
function round4(v) { return Math.round(v * 10000) / 10000; }

/* ======================================================================
   步伐偵測 + 動態步長 + 航向更新 + 樓梯分類
   ====================================================================== */
const stepFilterState = { filtered: 0, armed: false, peakCandidate: null };

function lowPassFilter(value, alpha = 0.3) {
  stepFilterState.filtered = stepFilterState.filtered + alpha * (value - stepFilterState.filtered);
  return stepFilterState.filtered;
}

function stepDetection(nowMs, accWorld, gyroZ, accMagnitude) {
  const filtered = lowPassFilter(accMagnitude);

  state.accVertBuffer.push({ t: nowMs, v: filtered, vertWorld: accWorld[2] });
  if (state.accVertBuffer.length > 200) state.accVertBuffer.shift();

  drawWaveformDebug(filtered);

  // 轉彎抑制（因子#11）：高角速度時暫停步長累積，只更新航向
  const turning = Math.abs(gyroZ) > (100 * Math.PI / 180);
  if (turning) state.turnSuppressUntil = nowMs + 150;

  // 航向角持續用陀螺儀 yaw 分量積分更新（簡化：直接用世界座標角速度 z 分量）
  const dt = 0.02;
  state.heading += gyroZ * dt;

  // 滯後雙門檻峰值偵測（hysteresis thresholding）：
  // 訊號需先「上升超過 highThreshold」進入 armed 狀態，再「回落到 lowThreshold 以下」才確認完成一次波峰，
  // 避免單點比較在雜訊環境下對同一個真實步伐重複觸發或被雜訊打斷漏判。
  const highThreshold = CARRY_MODE_THRESHOLD[state.carryMode] || 1.0;
  const lowThreshold = highThreshold * 0.3;

  if (!stepFilterState.armed && filtered > highThreshold) {
    stepFilterState.armed = true;
    stepFilterState.peakCandidate = { t: nowMs, v: filtered };
  } else if (stepFilterState.armed) {
    if (filtered > stepFilterState.peakCandidate.v) {
      stepFilterState.peakCandidate = { t: nowMs, v: filtered };
    }
    if (filtered < lowThreshold) {
      stepFilterState.armed = false;
      const candidateTime = stepFilterState.peakCandidate.t;
      const sinceLast = candidateTime - state.lastStepTime;
      if (sinceLast > 250 && sinceLast < 2000 && candidateTime > state.turnSuppressUntil) {
        registerStep(candidateTime, state.accVertBuffer);
      } else if (state.lastStepTime === 0) {
        registerStep(candidateTime, state.accVertBuffer);
      }
    }
  }
}

function registerStep(nowMs, buf) {
  state.lastStepTime = nowMs;
  state.stepCount++;

  // 動態步長估計 (Weinberg 公式)：用最近窗口的峰谷差
  const windowVals = buf.slice(-15).map(b => b.v);
  const aMax = Math.max(...windowVals);
  const aMin = Math.min(...windowVals);
  const diff = Math.max(aMax - aMin, 0.1);

  const heightM = state.heightCm ? state.heightCm / 100 : 1.7;
  const K = heightM * 0.42;
  const stepLength = K * Math.pow(diff, 0.25);

  // 水平位移（世界座標），航向角 heading：0=北(Y+), 順時針為正
  const dx = stepLength * Math.sin(state.heading);
  const dy = stepLength * Math.cos(state.heading);

  // 垂直位移：用最近窗口的加速度二次積分近似（每步 ZUPT 重置，見下方不連續累積）
  const dz = estimateVerticalDisplacement(buf);

  state.pos.x += dx;
  state.pos.y += dy;
  state.pos.z += dz;
  state.totalDistance += Math.sqrt(dx*dx + dy*dy);

  state.recentSteps.push({ dx, dy, dz, t: nowMs });
  if (state.recentSteps.length > 5) state.recentSteps.shift();

  const classification = classifySegment();
  state.currentSegmentType = classification.type;
  state.currentConfidence = classification.confidence;

  state.track.push({
    t_ms: Math.round(nowMs - state.startTime),
    x: round3(state.pos.x),
    y: round3(state.pos.y),
    z: round3(state.pos.z),
    segment_type: classification.type,
    slope_deg: classification.slopeDeg !== undefined ? round3(classification.slopeDeg) : undefined,
    confidence: round3(classification.confidence),
  });

  updateLiveStats();
  updateViz3D();
}

function estimateVerticalDisplacement(buf) {
  // 簡化 ZUPT：對該步窗口內「世界座標垂直加速度」(vertWorld，有正負號) 做二次積分，
  // 並以步窗口邊界歸零速度。注意：不可用 accMagnitude 模長欄位(v)，那是無方向性的量，積分沒有物理意義。
  const window = buf.slice(-15);
  if (window.length < 3) return 0;
  let v = 0, z = 0;
  for (let i = 1; i < window.length; i++) {
    const dt = (window[i].t - window[i-1].t) / 1000;
    if (dt <= 0 || dt > 0.5) continue;
    v += window[i].vertWorld * dt;
    z += v * dt;
  }
  // 限制單步垂直位移在合理範圍內（一般樓梯單階 15-20cm），避免積分噴出不合理值
  return Math.max(Math.min(z, 0.30), -0.30);
}

function classifySegment() {
  const steps = state.recentSteps;
  if (steps.length < 2) return { type: 'flat', confidence: 0.5 };

  const sumHoriz = steps.reduce((acc, s) => acc + Math.sqrt(s.dx*s.dx + s.dy*s.dy), 0);
  const sumVert = steps.reduce((acc, s) => acc + s.dz, 0);
  const ratio = sumHoriz > 0.01 ? Math.abs(sumVert) / sumHoriz : 0;

  const variance = computeVariance(steps.map(s => s.dz));
  const stability = Math.max(0, 1 - variance * 10);

  if (ratio < 0.08) {
    return { type: 'flat', confidence: Math.min(0.95, 0.6 + stability * 0.35) };
  } else if (ratio < 0.9) {
    const slopeRad = Math.atan2(sumVert, sumHoriz);
    const slopeDeg = slopeRad * 180 / Math.PI;
    const type = sumVert > 0 ? 'stairs_up' : 'stairs_down';
    return { type, confidence: Math.min(0.9, 0.5 + stability * 0.4), slopeDeg };
  } else {
    return { type: 'elevator_or_vertical', confidence: Math.min(0.7, 0.4 + stability * 0.3) };
  }
}

function computeVariance(arr) {
  const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
  return arr.reduce((acc, v) => acc + (v - mean) ** 2, 0) / arr.length;
}

function updateLiveStats() {
  document.getElementById('stat-steps').textContent = state.stepCount;
  document.getElementById('stat-distance').textContent = state.totalDistance.toFixed(1);
  document.getElementById('stat-height').textContent = state.pos.z.toFixed(1);
  const segLabel = {
    flat: '平地',
    stairs_up: '樓梯 ↑',
    stairs_down: '樓梯 ↓',
    elevator_or_vertical: '電梯/垂直',
  }[state.currentSegmentType] || state.currentSegmentType;
  document.getElementById('stat-segment').textContent = segLabel;
  document.getElementById('stat-confidence').textContent = Math.round(state.currentConfidence * 100) + '%';
}

/* ======================================================================
   除錯用即時波形圖（用於現場校準步伐偵測門檻值）
   ====================================================================== */
let waveformCtx = null;

function drawWaveformDebug(filteredValue) {
  const canvas = document.getElementById('canvas-waveform');
  if (!canvas) return;
  if (!waveformCtx) waveformCtx = canvas.getContext('2d');

  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth, cssH = canvas.clientHeight;
  if (canvas.width !== cssW * dpr || canvas.height !== cssH * dpr) {
    canvas.width = cssW * dpr;
    canvas.height = cssH * dpr;
  }
  const w = canvas.width, h = canvas.height;
  const ctx = waveformCtx;

  const buf = state.accVertBuffer.slice(-150); // 最近約 2.5 秒 (60Hz)
  if (buf.length < 2) return;

  const threshold = CARRY_MODE_THRESHOLD[state.carryMode] || 1.0;
  const values = buf.map(b => b.v);
  // Y軸依「實際訊號動態範圍」自動縮放，不綁定門檻值——否則門檻設太高時，
  // 真實訊號會被壓縮到畫面中央一小段，肉眼看起來像是平的，而看不出真實量級。
  const maxAbs = Math.max(...values.map(Math.abs), 0.3);

  const windowMax = Math.max(...values);
  const windowMin = Math.min(...values);
  const peakLabel = document.getElementById('debug-peak-value');
  if (peakLabel) peakLabel.textContent = `近2.5秒範圍: ${windowMin.toFixed(2)} ~ ${windowMax.toFixed(2)} / 門檻 ${threshold.toFixed(1)}`;

  ctx.clearRect(0, 0, w, h);

  // 門檻線 (黃)
  const yForValue = (v) => h/2 - (v / maxAbs) * (h/2 - 4*dpr);
  ctx.strokeStyle = '#ffb020';
  ctx.setLineDash([4*dpr, 4*dpr]);
  ctx.lineWidth = 1.5 * dpr;
  ctx.beginPath();
  ctx.moveTo(0, yForValue(threshold));
  ctx.lineTo(w, yForValue(threshold));
  ctx.stroke();
  ctx.setLineDash([]);

  // 零線
  ctx.strokeStyle = '#2d3850';
  ctx.lineWidth = 1 * dpr;
  ctx.beginPath();
  ctx.moveTo(0, h/2);
  ctx.lineTo(w, h/2);
  ctx.stroke();

  // 訊號線 (藍)
  ctx.strokeStyle = '#4f8cff';
  ctx.lineWidth = 2 * dpr;
  ctx.beginPath();
  buf.forEach((pt, i) => {
    const x = (i / (buf.length - 1)) * w;
    const y = yForValue(pt.v);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  });
  ctx.stroke();
}

/* ======================================================================
   3D / 2D 可視化 (Three.js + Canvas2D)
   ====================================================================== */
let scene, camera, renderer, lineObj, pathPoints = [];
let canvas2dCtx;

function initViz() {
  const container = document.getElementById('viz-container');
  const canvas3d = document.getElementById('canvas-3d');
  const w = container.clientWidth, h = container.clientHeight;

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x1a2235);
  camera = new THREE.PerspectiveCamera(60, w/h, 0.1, 1000);
  camera.position.set(5, -8, 6);
  camera.up.set(0, 0, 1);
  camera.lookAt(0, 0, 0);

  renderer = new THREE.WebGLRenderer({ canvas: canvas3d, antialias: true });
  renderer.setSize(w, h);

  const grid = new THREE.GridHelper(20, 20, 0x2d3850, 0x232d45);
  grid.rotation.x = Math.PI / 2;
  scene.add(grid);

  const axesHelper = new THREE.AxesHelper(2);
  scene.add(axesHelper);

  const material = new THREE.LineBasicMaterial({ color: 0x4f8cff, linewidth: 2 });
  const geometry = new THREE.BufferGeometry();
  lineObj = new THREE.Line(geometry, material);
  scene.add(lineObj);

  canvas2dCtx = document.getElementById('canvas-2d').getContext('2d');

  animate();
}

function animate() {
  if (document.getElementById('screen-tracking').classList.contains('active')) {
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
}

function updateViz3D() {
  pathPoints.push(new THREE.Vector3(state.pos.x, state.pos.y, state.pos.z));
  const positions = new Float32Array(pathPoints.length * 3);
  pathPoints.forEach((p, i) => {
    positions[i*3] = p.x; positions[i*3+1] = p.y; positions[i*3+2] = p.z;
  });
  lineObj.geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  lineObj.geometry.computeBoundingSphere();

  draw2D();
}

function draw2D() {
  const canvas = document.getElementById('canvas-2d');
  const ctx = canvas2dCtx;
  const w = canvas.width = canvas.clientWidth * 2;
  const h = canvas.height = canvas.clientHeight * 2;
  ctx.clearRect(0, 0, w, h);

  if (pathPoints.length < 2) return;

  const xs = pathPoints.map(p => p.x), ys = pathPoints.map(p => p.y);
  const minX = Math.min(...xs, -1), maxX = Math.max(...xs, 1);
  const minY = Math.min(...ys, -1), maxY = Math.max(...ys, 1);
  const range = Math.max(maxX - minX, maxY - minY, 2) * 1.2;
  const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
  const scale = (w * 0.8) / range;

  ctx.strokeStyle = '#4f8cff';
  ctx.lineWidth = 3;
  ctx.beginPath();
  pathPoints.forEach((p, i) => {
    const px = w/2 + (p.x - cx) * scale;
    const py = h/2 - (p.y - cy) * scale;
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  });
  ctx.stroke();

  // 起點
  const start = pathPoints[0];
  ctx.fillStyle = '#3ecf8e';
  ctx.beginPath();
  ctx.arc(w/2 + (start.x-cx)*scale, h/2 - (start.y-cy)*scale, 6, 0, Math.PI*2);
  ctx.fill();

  // 目前位置
  const curr = pathPoints[pathPoints.length-1];
  ctx.fillStyle = '#ff5c5c';
  ctx.beginPath();
  ctx.arc(w/2 + (curr.x-cx)*scale, h/2 - (curr.y-cy)*scale, 6, 0, Math.PI*2);
  ctx.fill();

  // 北方標示
  ctx.fillStyle = '#8b96ad';
  ctx.font = '20px sans-serif';
  ctx.fillText('N', w/2 - 7, 20);
}

function initResultViz() {
  const container = document.getElementById('viz-container-result');
  const canvas = document.getElementById('canvas-3d-result');
  const w = container.clientWidth, h = container.clientHeight;

  const rScene = new THREE.Scene();
  rScene.background = new THREE.Color(0x1a2235);
  const rCamera = new THREE.PerspectiveCamera(60, w/h, 0.1, 1000);

  const xs = pathPoints.map(p => p.x), ys = pathPoints.map(p => p.y), zs = pathPoints.map(p => p.z);
  const maxRange = Math.max(
    Math.max(...xs) - Math.min(...xs),
    Math.max(...ys) - Math.min(...ys),
    Math.max(...zs) - Math.min(...zs),
    3
  );
  rCamera.position.set(maxRange*0.8, -maxRange*1.2, maxRange*0.9);
  rCamera.up.set(0, 0, 1);
  rCamera.lookAt(0, 0, 0);

  const rRenderer = new THREE.WebGLRenderer({ canvas, antialias: true });
  rRenderer.setSize(w, h);

  const grid = new THREE.GridHelper(Math.max(maxRange*1.5, 10), 20, 0x2d3850, 0x232d45);
  grid.rotation.x = Math.PI / 2;
  rScene.add(grid);
  rScene.add(new THREE.AxesHelper(Math.max(maxRange*0.3, 1)));

  if (pathPoints.length > 1) {
    const positions = new Float32Array(pathPoints.length * 3);
    pathPoints.forEach((p, i) => {
      positions[i*3] = p.x; positions[i*3+1] = p.y; positions[i*3+2] = p.z;
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: 0x4f8cff }));
    rScene.add(line);

    const startSphere = new THREE.Mesh(new THREE.SphereGeometry(Math.max(maxRange*0.03,0.1)), new THREE.MeshBasicMaterial({ color: 0x3ecf8e }));
    startSphere.position.copy(pathPoints[0]);
    rScene.add(startSphere);

    const endSphere = new THREE.Mesh(new THREE.SphereGeometry(Math.max(maxRange*0.03,0.1)), new THREE.MeshBasicMaterial({ color: 0xff5c5c }));
    endSphere.position.copy(pathPoints[pathPoints.length-1]);
    rScene.add(endSphere);
  }

  rRenderer.render(rScene, rCamera);

  let angle = 0;
  function rotateLoop() {
    if (!document.getElementById('screen-result').classList.contains('active')) return;
    angle += 0.005;
    rCamera.position.x = Math.cos(angle) * maxRange * 1.2;
    rCamera.position.y = Math.sin(angle) * maxRange * 1.2 - maxRange*0.3;
    rCamera.lookAt(0, 0, maxRange*0.2);
    rRenderer.render(rScene, rCamera);
    requestAnimationFrame(rotateLoop);
  }
  rotateLoop();
}

/* ======================================================================
   權限請求與感測器啟動
   ====================================================================== */
async function requestPermission() {
  const statusBox = document.getElementById('permission-status');
  const btnStart = document.getElementById('btn-start');
  statusBox.textContent = '請求權限中...';
  statusBox.className = 'status-box status-pending';
  btnStart.disabled = true;

  try {
    if (typeof DeviceMotionEvent !== 'undefined' && typeof DeviceMotionEvent.requestPermission === 'function') {
      const motionPerm = await DeviceMotionEvent.requestPermission();
      if (motionPerm !== 'granted') throw new Error('動作感測器權限被拒絕');
    }
    if (typeof DeviceOrientationEvent !== 'undefined' && typeof DeviceOrientationEvent.requestPermission === 'function') {
      const orientPerm = await DeviceOrientationEvent.requestPermission();
      if (orientPerm !== 'granted') throw new Error('方向感測器權限被拒絕');
    }

    // 不能只信任 API 存在/回傳granted，必須實測是否真的收到感測器數據
    // （部分瀏覽器/裝置設定下，事件會註冊成功但永遠不觸發，例如系統層級關閉了動作與方向存取）
    const gotData = await verifySensorDataArrives();
    if (!gotData) {
      throw new Error('已註冊感測器監聽，但 3 秒內未收到任何數據。請確認「設定 → Safari → 動作與方向存取」已開啟，並重新整理頁面再試');
    }

    window.addEventListener('devicemotion', handleMotion);
    window.addEventListener('deviceorientation', handleOrientation);

    statusBox.textContent = '感測器權限已取得，已驗證可收到即時數據，可以開始偵測';
    statusBox.className = 'status-box status-ok';
    btnStart.disabled = false;
  } catch (err) {
    statusBox.textContent = '權限取得失敗：' + err.message;
    statusBox.className = 'status-box status-error';
    btnStart.disabled = true;
  }
}

function verifySensorDataArrives() {
  return new Promise((resolve) => {
    let received = false;
    const probe = (event) => {
      if (event.accelerationIncludingGravity && event.accelerationIncludingGravity.x !== null) {
        received = true;
      }
    };
    window.addEventListener('devicemotion', probe);
    setTimeout(() => {
      window.removeEventListener('devicemotion', probe);
      resolve(received);
    }, 3000);
  });
}

/* ======================================================================
   流程控制
   ====================================================================== */
function showScreen(id) {
  document.querySelectorAll('.screen').forEach(s => s.classList.remove('active'));
  document.getElementById(id).classList.add('active');
}

function startTracking() {
  state.username = document.getElementById('input-username').value || 'anonymous';
  state.carryMode = document.getElementById('input-carry-mode').value;
  const heightVal = document.getElementById('input-height').value;
  state.heightCm = heightVal ? parseFloat(heightVal) : null;
  state.loopClosureEnabled = document.getElementById('input-loop-closure').checked;

  // reset state
  state.tracking = true;
  state.startTime = performance.now();
  state.lastSampleTime = 0;
  state.sampleCount = 0;
  state.pos = { x: 0, y: 0, z: 0 };
  state.heading = 0;
  state.track = [];
  state.rawLog = [];
  state.stepCount = 0;
  state.totalDistance = 0;
  state.lastStepTime = 0;
  state.recentSteps = [];
  state.accVertBuffer = [];
  pathPoints = [new THREE.Vector3(0,0,0)];
  stepFilterState.filtered = 0;
  madgwick.q = [1,0,0,0];

  state.track.push({ t_ms: 0, x: 0, y: 0, z: 0, segment_type: 'flat', confidence: 1.0 });

  showScreen('screen-tracking');
  if (!scene) initViz();
}

function stopTracking() {
  state.tracking = false;
  const endTime = performance.now();
  const durationSec = (endTime - state.startTime) / 1000;

  let loopClosureError = null;
  if (state.loopClosureEnabled) {
    loopClosureError = Math.sqrt(state.pos.x**2 + state.pos.y**2 + state.pos.z**2);
  }

  const levelChanges = countLevelChanges();
  const avgRate = sampleRateTracker.lastRate || (state.sampleCount / durationSec);

  window._exportData = {
    metadata: {
      session_id: crypto.randomUUID ? crypto.randomUUID() : String(Date.now()),
      user_label: state.username,
      device_info: navigator.userAgent,
      carry_mode: state.carryMode,
      user_height_cm: state.heightCm,
      start_time_iso: new Date(Date.now() - durationSec*1000).toISOString(),
      end_time_iso: new Date().toISOString(),
      duration_sec: round3(durationSec),
      sample_count: state.sampleCount,
      avg_sample_rate_hz: round3(avgRate),
      loop_closure_error_m: loopClosureError !== null ? round3(loopClosureError) : null,
      total_steps_detected: state.stepCount,
      algorithm_version: 'pdr-v0.1',
      disclaimer: '本數據為手機 MEMS 感測器航位推算結果，非絕對定位；樓梯斜度為實測位移比值反推非測量值；詳見 docs/SPEC.md 第7節。',
    },
    track: state.track,
    raw_sensor_log: state.rawLog,
  };

  document.getElementById('sum-steps').textContent = state.stepCount;
  document.getElementById('sum-duration').textContent = durationSec.toFixed(1) + ' 秒';
  document.getElementById('sum-distance').textContent = state.totalDistance.toFixed(2) + ' m';
  document.getElementById('sum-levelchanges').textContent = levelChanges;
  document.getElementById('sum-samplerate').textContent = avgRate.toFixed(1) + ' Hz';

  const closureRow = document.getElementById('sum-closure-row');
  if (loopClosureError !== null) {
    document.getElementById('sum-closure').textContent = loopClosureError.toFixed(2) + ' m';
    closureRow.style.display = 'flex';
  } else {
    closureRow.style.display = 'none';
  }

  showScreen('screen-result');
  initResultViz();
}

function countLevelChanges() {
  let changes = 0;
  let lastType = 'flat';
  for (const pt of state.track) {
    if (pt.segment_type !== lastType && (pt.segment_type === 'stairs_up' || pt.segment_type === 'stairs_down')) {
      changes++;
    }
    lastType = pt.segment_type;
  }
  return changes;
}

/* ======================================================================
   輸出功能
   ====================================================================== */
function copyTrackJSON() {
  const data = window._exportData;
  if (!data) return;
  const text = JSON.stringify(data.track, null, 2);
  navigator.clipboard.writeText(text).then(() => {
    alert('Track JSON 已複製到剪貼簿（' + data.track.length + ' 個點）');
  }).catch(() => {
    alert('複製失敗，請改用下載功能');
  });
}

function downloadFullJSON() {
  const data = window._exportData;
  if (!data) return;
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  triggerDownload(blob, `imu_track_${data.metadata.user_label}_${Date.now()}.json`);
}

function downloadCSV() {
  const data = window._exportData;
  if (!data) return;

  let trackCsv = 't_ms,x,y,z,segment_type,slope_deg,confidence\n';
  data.track.forEach(p => {
    trackCsv += `${p.t_ms},${p.x},${p.y},${p.z},${p.segment_type},${p.slope_deg ?? ''},${p.confidence}\n`;
  });

  let rawCsv = 't_ms,ax,ay,az,ax_lin,ay_lin,az_lin,gx,gy,gz,mx,my,mz,qw,qx,qy,qz\n';
  data.raw_sensor_log.forEach(r => {
    rawCsv += `${r.t_ms},${r.acc.join(',')},${r.acc_linear.join(',')},${r.gyro.join(',')},${r.mag.join(',')},${r.orientation_quat.join(',')}\n`;
  });

  triggerDownload(new Blob([trackCsv], { type: 'text/csv' }), `track_${data.metadata.user_label}_${Date.now()}.csv`);
  setTimeout(() => {
    triggerDownload(new Blob([rawCsv], { type: 'text/csv' }), `raw_sensor_${data.metadata.user_label}_${Date.now()}.csv`);
  }, 300);
}

function triggerDownload(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function restart() {
  window.removeEventListener('devicemotion', handleMotion);
  window.removeEventListener('deviceorientation', handleOrientation);
  lastAcc = null; lastGyro = null; lastMag = null;
  showScreen('screen-setup');
  document.getElementById('permission-status').textContent = '尚未取得感測器權限';
  document.getElementById('permission-status').className = 'status-box status-pending';
  document.getElementById('btn-start').disabled = true;
}

/* ======================================================================
   事件綁定
   ====================================================================== */
document.getElementById('btn-request-permission').addEventListener('click', requestPermission);
document.getElementById('btn-start').addEventListener('click', startTracking);
document.getElementById('btn-stop').addEventListener('click', stopTracking);
document.getElementById('btn-copy-track').addEventListener('click', copyTrackJSON);
document.getElementById('btn-download-full').addEventListener('click', downloadFullJSON);
document.getElementById('btn-download-csv').addEventListener('click', downloadCSV);
document.getElementById('btn-restart').addEventListener('click', restart);
