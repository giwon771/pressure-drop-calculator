import streamlit as st
import json
import math
import pandas as pd

# --- 데이터 로드 함수 ---
def load_json(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    except:
        return None

def interpolate(temp, properties, key):
    # properties가 리스트인 경우 (현재 JSON 구조)
    try:
        # 온도 기준 정렬
        sorted_props = sorted(properties, key=lambda x: x['temp'])
        temps = [float(p['temp']) for p in sorted_props]
        
        # 범위 밖 처리
        if temp <= temps[0]: return sorted_props[0][key]
        if temp >= temps[-1]: return sorted_props[-1][key]
        
        # 선형 보간
        for i in range(len(temps)-1):
            t1, t2 = temps[i], temps[i+1]
            if t1 <= temp <= t2:
                v1 = sorted_props[i][key]
                v2 = sorted_props[i+1][key]
                return v1 + (v2 - v1) * (temp - t1) / (t2 - t1)
    except Exception as e:
        return 0.0
    return 0.0
    
# --- 경제적 최적 지름 계산 엔진 (예제 4.6 로직) ---
def solve_economic_diameter(rho, mu, m_dot, c1, c2, t, n, a, b, f_multiplier, eta, epsilon):
    # 초기 가정값: 예제 4.6처럼 임의의 지름(4cm)으로 시작 [cite: 298]
    D_guess = 0.04  
    tolerance = 0.00001 # 0.01mm 오차 이내 수렴 조건
    max_iter = 50
    
    for i in range(max_iter):
        # A. 현재 지름으로 유속 및 Re 계산 [cite: 293]
        v = (4 * m_dot) / (rho * math.pi * D_guess**2)
        re = (4 * m_dot) / (math.pi * D_guess * mu)
        
        # B. 마찰 계수 f 산출 (Haaland 식 - Moody 다이어그램 대용) [cite: 302]
        if re > 2300:
            term = (epsilon / D_guess / 3.7)**1.11 + (6.9 / re)
            f = (-1.8 * math.log10(term))**-2
        else:
            f = 64 / re
            
        # C. 식 (4.12)를 이용해 새로운 D_opt 계산 [cite: 212, 285]
        # 에너지 단가 c2 단위를 $/kWh에서 W단위로 보정 [cite: 262]
        numerator = 40 * f * (m_dot**3) * (c2 / 1000) * t 
        denominator = n * (a + b) * (1 + f_multiplier) * c1 * eta * (math.pi**2) * (rho**2)
        
        D_new = (numerator / denominator)**(1 / (n + 5))
        
        # D. 수렴 확인: 가정한 D와 계산된 D의 차이가 작으면 종료 [cite: 312]
        if abs(D_new - D_guess) < tolerance:
            return D_new, f, re
        
        D_guess = D_new # 다음 루프를 위해 업데이트
        
    return D_guess, f, re

# --- 페이지 설정 ---
st.set_page_config(page_title="공학용 유체 설계 시스템 v8.0", layout="wide")

st.title("🚀 공학용 유체 수송 설계 시스템 (Web v8.0)")
st.markdown("### 배관 내 유동 분석 및 압력 손실 계산기")

# 데이터 로드
f_db = load_json('fluids_db.json')
p_db = load_json('pipe_db.json')

if not f_db or not p_db:
    st.error("데이터 파일(JSON)을 찾을 수 없습니다.")
    st.stop()

# --- 사이드바: 유체 및 온도 설정 [1] ---
with st.sidebar:
    st.header("[1] 유체 물성 설정")
    fluid_name = st.selectbox("대상 유체 선택", [f['name'] for f in f_db['fluids']])
    
    fix_temp = st.checkbox("상온 고정 (20°C)", value=True)
    if fix_temp:
        target_temp = 20.0
    else:
        target_temp = st.number_input("운전 온도 (°C)", value=20.0)

    fluid_data = next(item for item in f_db['fluids'] if item['name'] == fluid_name)
    rho = interpolate(target_temp, fluid_data['properties'], 'rho')
    mu = interpolate(target_temp, fluid_data['properties'], 'mu')
    
    st.info(f"선택 유체: {fluid_name}\n\n밀도: {rho:.2f} kg/m³\n\n점도: {mu:.6f} Pa·s")

# --- 메인 화면: 입력창 구성 [2] & [3] ---
col1, col2 = st.columns(2)

with col1:
    st.subheader("[2] 배관 규격 및 수치")
    
    # 배관 규격 선택
    nps_list = [p['nps'] for p in p_db['pipe_standards']]
    sel_nps = st.selectbox("NPS 선택", nps_list)
    sel_sch = st.selectbox("Schedule 선택", ["Sch 40", "Sch 80", "Sch 160"])
    
    # 자동 안지름 계산 (선택 시 기본값 제안)
    # 실제 구현 시에는 p_db에서 매칭되는 ID를 가져오는 로직 필요
    d_val = st.number_input("관 안지름(ID)", value=52.5, format="%.2f")
    d_unit = st.selectbox("직경 단위", ["mm", "m", "inch"], index=0)
    
    # 거리 입력
    l_val = st.number_input("배관 직선 거리(L)", value=10.0)
    l_unit = st.selectbox("거리 단위", ["m", "km"], index=0)

with col2:
    st.subheader("[3] 유동 및 부속품")
    
    # 유량/유속 입력
    v_val = st.number_input("유속 또는 유량 입력", value=1.0)
    v_unit = st.selectbox("입력 단위 선택", ["m/s", "m³/s", "L/min", "L/s", "BTU/hr (Water)"])
    
    st.markdown("---")
    # 부속품
    n_elbow = st.number_input("엘보(90°) 개수", min_value=0, value=0)
    n_valve = st.number_input("게이트밸브 개수", min_value=0, value=0)

# --- 계산 실행 버튼 ---
if st.button("🚀 설계 분석 및 계산 실행", use_container_width=True):
    # [1. 단위 변환]
    D = d_val / 1000 if d_unit == "mm" else (d_val * 0.0254 if d_unit == "inch" else d_val)
    L = l_val * 1000 if l_unit == "km" else l_val
    area = math.pi * (D**2) / 4

    # [2. 유속 변환]
    v_calc_process = ""
    if v_unit == "m/s":
        v = v_val
    elif v_unit == "m³/s":
        v = v_val / area
    elif v_unit == "L/min":
        v = (v_val / 60000) / area
    elif v_unit == "BTU/hr (Water)":
        q_m3s = (v_val / 500) * 0.00006309
        v = q_m3s / area
        v_calc_process = "(BTU/hr 기반 환산)"

    # [3. Re 및 흐름 판정]
    re = (rho * v * D) / mu
    
    if re <= 2300:
        flow_type, f, color = "층류 (Laminar)", 64/re, "blue"
    elif 2300 < re < 4000:
        flow_type, f, color = "⚠️ 천이 영역 (Transition)", 0.316*(re**-0.25), "orange"
        st.warning(f"경고: Reynolds Number({re:.1f})가 천이 영역에 있습니다. 유동이 불안정합니다.")
    else:
        flow_type, f, color = "난류 (Turbulent)", 0.316*(re**-0.25), "red"

    # [4. 압력 손실 및 BTU 환산]
    sum_k = (n_elbow * 0.9) + (n_valve * 0.2)
    dp_major = f * (L / D) * (rho * v**2 / 2)
    dp_minor = sum_k * (rho * v**2 / 2)
    dp_total = dp_major + dp_minor
    
    energy_loss_btu_lb = (dp_total / rho) * 0.00042906

    # [새로운 경제성 분석 호출 부분]
    # 사이드바에서 입력받은 변수들을 함수에 전달합니다.
    d_opt_m, final_f, final_re = solve_economic_diameter(
        rho, mu, m_dot, c1_value, c2, t_year, n_exponent, ann_a, ann_b, cost_f, eff_pump, 0.000046
    )
    
    d_opt_mm = d_opt_m * 1000
    
    # --- 결과 출력 ---
    st.subheader("💰 최적 경제적 직경 분석 결과")
    st.success(f"이 시스템의 최적 경제적 지름(D_opt)은 **{d_opt_mm:.2f} mm** 입니다.")
    st.info(f"선정 근거: 연간 운영비(에너지)와 초기 투자비의 합산 최소화 [cite: 14, 203]")
    res_col1, res_col2 = st.columns(2)
    
    with res_col1:
        st.metric("Reynolds Number (Re)", f"{re:.1f}")
        st.metric("흐름 상태", flow_type)
        
    with res_col2:
        st.metric("총 압력 손실 (ΔP)", f"{dp_total/1000:.4f} kPa")
        st.metric("에너지 손실", f"{energy_loss_btu_lb:.6f} BTU/lb")

    with st.expander("📝 상세 계산 리포트 보기"):
        st.code(f"""
[상세 설계 리포트]
- 유체: {fluid_name} ({target_temp} °C)
- 관경(D): {D:.4f} m / 유속(v): {v:.4f} m/s
- 마찰계수(f): {f:.4f} ({flow_type})
- 주손실: {dp_major:.2f} Pa / 부손실: {dp_minor:.2f} Pa
- 최종 ΔP: {dp_total:.2f} Pa
        """)
