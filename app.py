# --- [계산 및 결과 섹션: 기존 코드의 167번 줄부터 이 내용으로 덮어쓰세요] ---

if st.button("🚀 설계 분석 및 계산 실행", use_container_width=True):
    # 단위 환산 로직
    D_input = d_val / 1000 if d_unit == "mm" else (d_val * 0.0254 if d_unit == "inch" else d_val)
    L = l_val * 1000 if l_unit == "km" else l_val
    area = math.pi * (D_input**2) / 4

    if v_unit == "m/s": v = v_val
    elif v_unit == "m³/s": v = v_val / area
    elif v_unit == "L/min": v = (v_val / 60000) / area
    else: v = (v_val / 1000) / area
    
    m_dot = rho * v * area 

    # 1. 이론적 최적 지름 계산 (시행착오법)
    d_opt_m, final_f, final_re = solve_economic_diameter(
        rho, mu, m_dot, c1_value, c2, t_year, n_exponent, ann_a, ann_b, cost_f, eff_pump, 0.000046
    )

    # 2. 상용 배관 규격 추천 로직 (Database Matching)
    recommended_pipe = None
    min_diff = float('inf')
    for p in p_db['pipe_standards']:
        if sel_sch in p['schedules']:
            db_id = p['schedules'][sel_sch]['id'] / 1000 # m 환산
            diff = abs(db_id - d_opt_m)
            if diff < min_diff:
                min_diff = diff
                recommended_pipe = {"nps": p['nps'], "id": db_id}

    # 3. 추천 배관 기반 경제성 재계산
    D_real = recommended_pipe['id']
    v_real = (4 * m_dot) / (rho * math.pi * D_real**2)
    re_real = (rho * v_real * D_real) / mu
    
    # 실제 마찰계수 산출 (Haaland 식 적용)
    term = (0.000046 / D_real / 3.7)**1.11 + (6.9 / re_real)
    f_real = (-1.8 * math.log10(term))**-2 if re_real > 2300 else 64/re_real
    
    # 실제 구간 압력 손실 및 연간 운영 비용 (OPEX)
    delta_p = f_real * (L / D_real) * (rho * v_real**2 / 2)
    power_kw = (delta_p * (m_dot / rho)) / (eff_pump * 1000)
    annual_energy_cost = power_kw * t_year * c2
    
    # 연간 고정 비용 (CAPEX 환산)
    annual_fixed_cost = c1_value * (D_real**n_exponent) * (ann_a + ann_b) * (1 + cost_f) * L
    total_annual_cost = annual_energy_cost + annual_fixed_cost

    st.divider()

    # --- [섹션 1: 이론적 D_opt 도출 과정 (기존 선호 화면)] ---
    st.subheader("💰 1. 경제적 최적 지름($D_{opt}$) 도출 결과")
    st.success(f"이론적 수치 해석 결과, 최적 경제 지름은 **{d_opt_m*1000:.2f} mm** 입니다.")
    
    with st.expander("🔍 수치 해석 및 시행착오법 수렴 리포트", expanded=True):
        st.markdown(f"""
        ### 최적화 알고리즘 (Darby's Equation)
        본 시스템은 **시행착오법(Trial and Error)**을 통해 아래 비선형 방정식을 만족하는 지점을 수렴시켰습니다.
        """)
        # LaTeX 수식 출력
        st.latex(rf"D_{{opt}} = \left[ \frac{{40 \cdot {final_f:.4f} \cdot {m_dot:.2f}^3 \cdot {c2/1000:.6f} \cdot {t_year}}}{{{n_exponent} \cdot ({ann_a:.3f} + {ann_b:.2f}) \cdot (1 + {cost_f:.1f}) \cdot {c1_value} \cdot {eff_pump} \cdot \pi^2 \cdot {rho:.0f}^2}} \right]^{{\frac{{1}}{{{n_exponent}+5}}}}")
        
        st.markdown(f"""
        **최종 수렴 데이터:**
        - **수렴 마찰 계수 ($f$):** {final_f:.5f}
        - **설계 Reynolds 수 ($Re$):** {final_re:.1f}
        """)

    # --- [섹션 2: 상용 규격 추천 및 비용 분석 (신규 피드백 반영)] ---
    st.subheader("📋 2. 상용 규격 추천 및 연간 총 비용(TAC)")
    
    c_res1, c_res2 = st.columns(2)
    with c_res1:
        st.info(f"**추천 규격:** NPS {recommended_pipe['nps']} (Sch.{sel_sch})")
        st.write(f"- 실제 배관 안지름: {D_real*1000:.2f} mm")
        st.write(f"- 이론치와의 편차: {min_diff*1000:.4f} mm")
        
    with c_res2:
        st.metric("총 연간 비용 (TAC)", f"$ {total_annual_cost:,.2f} /yr")
        st.caption(f"운전비(OPEX): ${annual_energy_cost:,.0f} | 시설비(CAPEX): ${annual_fixed_cost:,.0f}")

    # 비용 구성 시각화
    st.write("### 연간 소요 비용 구성 (Annual Cost Breakdown)")
    cost_data = {"에너지 운영 비용": annual_energy_cost, "시설 설치비(연간화)": annual_fixed_cost}
    st.bar_chart(cost_data)

    # Re 판정 결과 메트릭 (하단 배치)
    c1_res, c2_res = st.columns(2)
    
    # 현재 유동 상태 판정
    re_current = (rho * v_real * D_real) / mu
    if re_current <= 2300: flow_type = "층류 (Laminar)"
    elif 2300 < re_current < 4000: flow_type = "천이 영역 (Transitional)"
    else: flow_type = "난류 (Turbulent)"
    
    c1_res.metric("설계 Re 수 (실제 배관 기준)", f"{re_current:.1f}")
    c2_res.metric("현재 유동 상태", flow_type)
