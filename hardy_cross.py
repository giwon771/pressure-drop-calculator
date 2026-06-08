import streamlit as st
import math
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# --- 1. 파이프 객체 정의 ---
class Pipe:
    def __init__(self, id, node_start, node_end, length, diameter, roughness, initial_q=0.0):
        self.id = id
        self.start = node_start
        self.end = node_end
        self.L = length
        self.D = diameter
        self.epsilon = roughness
        self.Q = initial_q  
        self.rho = 1000
        self.mu = 0.89e-3
        
    def get_c_constant(self):
        return (8 * self.rho * self.L) / (math.pi**2 * (self.D**5))

    def get_friction_factor(self):
        re = (4 * self.rho * abs(self.Q)) / (math.pi * self.D * self.mu)
        if re < 2300: 
            return 64 / re if re > 0 else 0.02
        term1 = self.epsilon / (3.7065 * self.D)
        term2 = 5.0452 / re
        term3 = (1 / 2.8257) * (self.epsilon / self.D)**1.1098 + 5.8506 / (re**0.8981)
        return (-2.0 * math.log10(term1 - term2 * math.log10(term3)))**-2

    def get_delta_p(self):
        return self.get_c_constant() * self.get_friction_factor() * self.Q * abs(self.Q)

    def get_dp_over_q(self):
        return abs(self.get_delta_p() / self.Q) if self.Q != 0 else 0

    def get_delta_p_per_km(self):
        total_pa = abs(self.get_delta_p())
        if self.L == 0: return 0
        pa_per_m = total_pa / self.L
        kpa_per_km = (pa_per_m * 1000) / 1000  
        return kpa_per_km

# --- 2. 동적 하디 크로스 수치 해석 엔진 ---
def run_dynamic_hardy_cross(pipes, loops_nodes, max_iter=30, tolerance=1e-5):
    iteration_history = []
    pipe_dict = {(p.start, p.end): p for p in pipes}
    pipe_dict.update({(p.end, p.start): p for p in pipes})
    
    for idx in range(max_iter):
        dq_list = []
        for loop in loops_nodes:
            sum_dp = 0.0
            sum_dp_q = 0.0
            for i in range(len(loop)):
                n1 = loop[i]
                n2 = loop[(i + 1) % len(loop)]
                if (n1, n2) in pipe_dict:
                    p = pipe_dict[(n1, n2)]
                    sign = 1.0 if p.start == n1 else -1.0
                    sum_dp += sign * p.get_delta_p()
                    sum_dp_q += p.get_dp_over_q()
            dq = -sum_dp / (2 * sum_dp_q) if sum_dp_q != 0 else 0
            dq_list.append(dq)
            
        q_adjustments = {p.id: 0.0 for p in pipes}
        for l_idx, loop in enumerate(loops_nodes):
            dq = dq_list[l_idx]
            for i in range(len(loop)):
                n1 = loop[i]
                n2 = loop[(i + 1) % len(loop)]
                p = pipe_dict[(n1, n2)]
                sign = 1.0 if p.start == n1 else -1.0
                q_adjustments[p.id] += sign * dq
                
        for p in pipes:
            p.Q += q_adjustments[p.id]
            
        iteration_history.append(dq_list)
        if all(abs(dq) < tolerance for dq in dq_list):
            break
    return iteration_history

# --- 3. Streamlit 실행 메인 함수 ---
def run_hardy_cross():
    if 'pipe_data' not in st.session_state:
        st.session_state.pipe_data = [
            {"id": 1, "start": "A", "end": "B", "L": 300.0, "D": 0.250, "init_q": 0.060, "sx": 300.0, "sy": 0.0, "ex": 0.0, "ey": 0.0},
            {"id": 2, "start": "B", "end": "E", "L": 250.0, "D": 0.200, "init_q": 0.020, "sx": 0.0, "sy": 0.0, "ex": 0.0, "ey": 250.0},
            {"id": 3, "start": "E", "end": "D", "L": 300.0, "D": 0.200, "init_q": -0.040, "sx": 0.0, "sy": 250.0, "ex": 300.0, "ey": 250.0},
            {"id": 4, "start": "D", "end": "A", "L": 250.0, "D": 0.250, "init_q": -0.065, "sx": 300.0, "sy": 250.0, "ex": 300.0, "ey": 0.0},
            {"id": 5, "start": "B", "end": "C", "L": 300.0, "D": 0.200, "init_q": 0.040, "sx": 0.0, "sy": 0.0, "ex": -300.0, "ey": 0.0},
            {"id": 6, "start": "C", "end": "F", "L": 250.0, "D": 0.200, "init_q": 0.028, "sx": -300.0, "sy": 0.0, "ex": -300.0, "ey": 250.0},
            {"id": 7, "start": "F", "end": "E", "L": 300.0, "D": 0.150, "init_q": -0.035, "sx": -300.0, "sy": 250.0, "ex": 0.0, "ey": 250.0}
        ]

    temp_pos = {}
    for p in st.session_state.pipe_data:
        temp_pos[p['start']] = (p.get('sx', 0.0), p.get('sy', 0.0))
        temp_pos[p['end']] = (p.get('ex', 0.0), p.get('ey', 0.0))
    existing_nodes = sorted(list(temp_pos.keys()))

    st.sidebar.header("[1] 시스템 가동 조건 설정")
    if len(existing_nodes) > 0:
        control_node = st.sidebar.selectbox("펌프/밸브 제어 노드 선택", existing_nodes)
    else:
        st.sidebar.warning("⚠️ 등록된 노드가 없습니다.")
        control_node = None
        
    total_inflow = st.sidebar.slider("시스템 총 유입 유량 (m³/s)", 0.05, 0.40, 0.125, step=0.005)
    roughness_val = st.sidebar.number_input("관 절대 조도 (m)", value=0.00025, format="%.5f")
    pump_eff = st.sidebar.slider("펌프 효율 (η)", 0.5, 0.9, 0.75)

    st.subheader("🕹️ 스마트 배관망 그래픽 배치 조립 판넬")
    
    with st.expander("📐 좌표 입력 없이 방향 선택으로 배관 쉽게 연장하기", expanded=False):
        if len(existing_nodes) == 0:
            st.info("💡 현재 네트워크가 비어 있습니다. 첫 배관의 시작점(원점)을 배치합니다.")
            c1, c2, c3, c4 = st.columns(4)
            p_start = c1.text_input("시작 노드 이름", value="A").strip().upper()
            p_end = c2.text_input("끝 노드 이름", value="B").strip().upper()
            p_L = c3.number_input("실제 배관 길이 L (m)", min_value=1.0, value=200.0, key="blank_L")
            p_D = c4.number_input("직경 D (m)", min_value=0.01, value=0.200, key="blank_D")
            
            if st.button("🚀 최초 원점 배관 생성", use_container_width=True):
                st.session_state.pipe_data.append({
                    "id": 1, "start": p_start, "end": p_end, "L": p_L, "D": p_D, "init_q": 0.05,
                    "sx": 0.0, "sy": 0.0, "ex": p_L, "ey": 0.0
                })
                st.rerun()
        else:
            col_ui1, col_ui2 = st.columns(2)
            
            with col_ui1:
                st.markdown("**1. 출발점 설정**")
                p_start = st.selectbox("어느 노드에서 배관을 연장할까요?", existing_nodes)
                sx, sy = temp_pos[p_start]
                
                st.markdown("**2. 연장할 방향 모드 선택**")
                layout_mode = st.radio("가설 형태", ["직각 방향 가설", "대각선 각도 지정 가설"])
                
                if layout_mode == "직각 방향 가설":
                    direction = st.radio(
                        "어느 방향으로 관을 가설합니까?",
                        ["우측으로 연장 (+X)", "좌측으로 연장 (-X)", "위로 연장 (+Y)", "아래로 연장 (-Y)"]
                    )
                    angle = 0.0
                else:
                    angle = st.slider("배관 가설 각도 입력 (도, °)", -180, 180, 45, step=5)
                    fig_circle, ax_c = plt.subplots(figsize=(2.2, 2.2))
                    circle = plt.Circle((0,0), 1.0, color='#BDC3C7', fill=False, linestyle='--', linewidth=1.2)
                    ax_c.add_patch(circle)
                    ax_c.axhline(0, color='#BDC3C7', linewidth=0.8, linestyle=':')
                    ax_c.axvline(0, color='#BDC3C7', linewidth=0.8, linestyle=':')
                    rad_preview = math.radians(angle)
                    vx, vy = math.cos(rad_preview), math.sin(rad_preview)
                    ax_c.quiver(0, 0, vx, vy, angles='xy', scale_units='xy', scale=1, color='#E74C3C', width=0.07)
                    ax_c.text(vx*1.3, vy*1.3, f"{angle}°", color='#E74C3C', fontsize=9, weight='bold', ha='center', va='center')
                    ax_c.set_xlim(-1.5, 1.5)
                    ax_c.set_ylim(-1.5, 1.5)
                    ax_c.axis('off')
                    plt.tight_layout()
                    st.pyplot(fig_circle)
                
            with col_ui2:
                st.markdown("**3. 도달점 및 스펙 지정**")
                p_end = st.text_input("도달 노드 이름 입력 (이미 있는 노드면 자동 스냅)", value="G").strip().upper()
                p_L = st.number_input("배관 물리 길이 L (m)", min_value=1.0, value=200.0)
                p_D = st.number_input("배관 직경 D (m)", min_value=0.01, value=0.200, format="%.3f")
                p_q = st.number_input("초기 가정 유량 (m³/s)", value=0.020, format="%.3f")

            if layout_mode == "직각 방향 가설":
                if direction == "우측으로 연장 (+X)": ex, ey = sx + p_L, sy
                elif direction == "좌측으로 연장 (-X)": ex, ey = sx - p_L, sy
                elif direction == "위로 연장 (+Y)": ex, ey = sx, sy + p_L
                else: ex, ey = sx, sy - p_L
            else:
                rad = math.radians(angle)
                ex = sx + p_L * math.cos(rad)
                ey = sy + p_L * math.sin(rad)
            
            is_snap = False
            if p_end in temp_pos:
                ex, ey = temp_pos[p_end]
                is_snap = True

            if st.button("🛠️ 지정 방향으로 배관라인 즉시 가설", use_container_width=True):
                if p_start == p_end:
                    st.error("시작 노드와 끝 노드가 같으면 루프 연산이 불가합니다.")
                else:
                    new_id = len(st.session_state.pipe_data) + 1
                    st.session_state.pipe_data.append({
                        "id": new_id, "start": p_start, "end": p_end, "L": p_L, "D": p_D, "init_q": p_q,
                        "sx": sx, "sy": sy, "ex": ex, "ey": ey
                    })
                    st.rerun()

    df_pipes = pd.DataFrame(st.session_state.pipe_data)
    if not df_pipes.empty:
        st.markdown("📝 **배관 스펙 실시간 편집 테이블** (셀을 더블클릭하여 수치를 즉시 변경해 보세요)")
        edited_df = st.data_editor(
            df_pipes[["id", "start", "end", "L", "D", "init_q"]],
            disabled=["id", "start", "end"], 
            use_container_width=True,
            key="pipe_editor"
        )
        
        for index, row in edited_df.iterrows():
            st.session_state.pipe_data[index]["L"] = float(row["L"])
            st.session_state.pipe_data[index]["D"] = float(row["D"])
            st.session_state.pipe_data[index]["init_q"] = float(row["init_q"])
            
            p_item = st.session_state.pipe_data[index]
            if p_item["sx"] == p_item["ex"]: 
                if p_item["ey"] >= p_item["sy"]: p_item["ey"] = p_item["sy"] + p_item["L"]
                else: p_item["ey"] = p_item["sy"] - p_item["L"]
            elif p_item["sy"] == p_item["ey"]: 
                if p_item["ex"] >= p_item["sx"]: p_item["ex"] = p_item["sx"] + p_item["L"]
                else: p_item["ex"] = p_item["sx"] - p_item["L"]

        # --- 🗑️ [신규 기능 통합] 지정 배관 라인 타겟 개별 삭제 패널 ---
        with st.expander("❌ 불필요한 특정 배관라인 선택하여 제거하기", expanded=False):
            # 현재 살아있는 파이프 ID 리스트 추출
            pipe_ids = [p["id"] for p in st.session_state.pipe_data]
            pipe_to_delete = st.selectbox("철거할 파이프 번호 선택", pipe_ids, format_func=lambda x: f"Pipe #{x}")
            
            if st.button("🚨 선택한 배관라인 즉시 철거", use_container_width=True):
                # 선택한 ID 항목만 쏙 빼고 데이터 재구성
                st.session_state.pipe_data = [p for p in st.session_state.pipe_data if p["id"] != pipe_to_delete]
                
                # ID 순번 꼬임 방지를 위한 1번부터 순차 재정렬 작업(Index Re-ordering)
                for idx, p in enumerate(st.session_state.pipe_data):
                    p["id"] = idx + 1
                    
                st.toast(f"Pipe #{pipe_to_delete} 선로가 안전하게 철거되었습니다.")
                st.rerun()
    
    reset_col1, reset_col2 = st.columns(2)
    if reset_col1.button("🔄 교재 예제 데이터로 초기화", use_container_width=True):
        if 'pipe_data' in st.session_state: del st.session_state.pipe_data
        st.rerun()
    if reset_col2.button("🗑️ 전체 노드 삭제 (Blank Reset)", use_container_width=True):
        st.session_state.pipe_data = []
        st.rerun()

    if len(st.session_state.pipe_data) == 0:
        st.warning("⚠️ 현재 시스템 내에 배관 선로가 전혀 없습니다. 위쪽 판넬에서 최초 원점 배관을 생성해 주세요.")
        return

    G_setup = nx.Graph()
    pipes = []
    pos = {}  
    
    scale = total_inflow / 0.125 if len(st.session_state.pipe_data) == 7 else 1.0
    for p in st.session_state.pipe_data:
        pipes.append(Pipe(p['id'], p['start'], p['end'], p['L'], p['D'], roughness_val, p['init_q'] * scale))
        G_setup.add_edge(p['start'], p['end'])
        pos[p['start']] = (p.get('sx', 0.0), p.get('sy', 0.0))
        pos[p['end']] = (p.get('ex', 0.0), p.get('ey', 0.0))

    auto_loops = nx.cycle_basis(G_setup)
    st.info(f"🔍 시스템 위상 분석 완료: 감지된 독립 루프 개수 = **{len(auto_loops)}개**")

    if len(auto_loops) == 0:
        st.warning("⚠️ 현재 폐회로(Loop)가 완전히 닫히지 않은 트리형 배관입니다. 하디 크로스 해석을 수행하려면 노드와 노드를 이어 루프를 완성해 주세요.")
        fig, ax = plt.subplots(figsize=(7, 4.5))
        nx.draw_networkx_nodes(G_setup, pos, node_size=600, node_color='#EAEDED', ax=ax)
        nx.draw_networkx_labels(G_setup, pos, font_size=11, font_weight='bold', ax=ax)
        nx.draw_networkx_edges(G_setup, pos, width=2, edge_color='#7F8C8D', ax=ax)
        
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.set_xlabel("X-Axis Distance (m)", fontsize=10, fontweight='bold')
        ax.set_ylabel("Y-Axis Distance (m)", fontsize=10, fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6)
        st.pyplot(fig)
        return

    history = run_dynamic_hardy_cross(pipes, auto_loops)

    # --- UI 레이아웃 화면 표시 ---
    col1, col2 = st.columns([3, 2])
    
    with col1:
        st.subheader("🖼️ 실제 길이 스케일 반영 2D 배관 플랜 도면")
        G_draw = nx.DiGraph()
        for p in pipes:
            if p.Q >= 0: G_draw.add_edge(p.start, p.end, weight=abs(p.Q), id=p.id)
            else: G_draw.add_edge(p.end, p.start, weight=abs(p.Q), id=p.id)
            
        fig, ax = plt.subplots(figsize=(7, 4.5))
        nx.draw_networkx_nodes(G_draw, pos, node_size=600, node_color='#D6EAF8', ax=ax)
        nx.draw_networkx_labels(G_draw, pos, font_size=11, font_weight='bold', ax=ax)
        
        if control_node and (control_node in pos):
            nx.draw_networkx_nodes(G_draw, pos, nodelist=[control_node], node_size=700, node_color='#FF5733', ax=ax)
            
        edges = G_draw.edges(data=True)
        weights = [max(e[2]['weight'] * 150, 1.5) for e in edges]
        nx.draw_networkx_edges(G_draw, pos, width=weights, edge_color='#2C3E50', arrowsize=18, ax=ax)
        
        edge_labels = {}
        for p in pipes:
            k = (p.start, p.end) if p.Q >= 0 else (p.end, p.start)
            edge_labels[k] = f"#{p.id} ({p.L:.0f}m)\n{abs(p.Q):.3f}m³/s"
            
        nx.draw_networkx_edge_labels(G_draw, pos, edge_labels=edge_labels, font_size=8, font_color='red')
        
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.set_xlabel("X-Axis Distance (m)", fontsize=10, fontweight='bold')
        ax.set_ylabel("Y-Axis Distance (m)", fontsize=10, fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6)
        
        plt.tight_layout()
        st.pyplot(fig)

    with col2:
        st.subheader("📊 계통 마찰 손실 및 소요 동력 리포트")
        
        total_dp_loss = sum(abs(p.get_delta_p()) for p in pipes)
        required_power_w = total_dp_loss * total_inflow
        required_power_kw = required_power_w / (pump_eff * 1000)
        
        st.metric("계통 총 마찰 압력 강하 (Total ΔP)", f"{total_dp_loss:,.1f} N/m²")
        st.metric("이론적 소요 동력 (Required Power)", f"{required_power_kw:.2f} kW")
        
        st.write("---")
        st.markdown("### 🔌 Grundfos 지능형 수력 범위(Hydraulic Coverage) 추천 모듈")
        
        head_loss_ratio = total_dp_loss / total_inflow  
        
        if head_loss_ratio >= 400000 and total_inflow < 0.15:
            series_name = "CR 시리즈 (Vertical Multistage)"
            series_desc = "본 계통은 유량 대비 **마찰 손실 압력 강하가 매우 높은 고양정 환경**입니다. 따라서 컴팩트한 바닥 면적에서 임펠러를 수직 다단으로 배열하여 초고압 토출을 구현하는 그룬포스 **CR 시리즈**가 유체 기계학적으로 가장 최적입니다."
            pumps_list = [
                {"search_name": "CR 15-3", "model": "Grundfos CR 15-3 A-A-A-E-HQQE (정격 고압 수직 다단형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"},
                {"search_name": "CR 45-2", "model": "Grundfos CR 45-2 A-F-A-E-HQQE (대유량 커버 수직 다단형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"}
            ]
        elif total_inflow >= 0.25:
            series_name = "NK 시리즈 (Long-Coupled Volute)"
            series_desc = "본 계통은 독립 루프망 전체를 관통하는 **순환 유량이 대규모인 계통**입니다. 따라서 모터와 펌프 축이 커플링으로 분리되어 연속 대유량 공정에서 기계적 진동·소음을 원천 차단하는 그룬포스 **NK 대형 볼류트 시리즈**를 강력하게 제안합니다."
            pumps_list = [
                {"search_name": "NK 100-200", "model": "Grundfos NK 100-200/219 (산업용 대유량 장축 볼류트형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"},
                {"search_name": "NK 150-315", "model": "Grundfos NK 150-315/304 (플랜트 메인 대용량 볼류트형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"}
            ]
        else:
            series_name = "NB 시리즈 (End-Suction Standard)"
            series_desc = "본 계통은 유량과 압력 강하비가 **가장 경제적인 평형 균형을 이루는 범용 계통**입니다. 따라서 전 세계 플랜트 표준으로 가장 널리 쓰이며 흡입구와 토출구가 직각을 이루어 유지보수 비용이 가장 저렴한 **NB 엔드석션 시리즈**가 최적입니다."
            pumps_list = [
                {"search_name": "NB 50-160", "model": "Grundfos NB 50-160/154 (표준 단단 엔드석션형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"},
                {"search_name": "NB 80-160", "model": "Grundfos NB 80-160/177 (대유량 고속 엔드석션형)", "spec": f"추천 정격 동력: {required_power_kw*1.15:.1f} kW 내외 | 60Hz 3상"}
            ]

        st.info(f"🧭 **계통 분석 진단 결과: {series_name} 매칭**")
        st.caption(series_desc)
        
        for idx, pump in enumerate(pumps_list):
            unique_key = f"p_map_btn_{required_power_kw:.2f}_{idx}_{pump['search_name']}"
            with st.container(border=True):
                st.markdown(f"**🏅 권장 세부 계열 기종 #{idx+1}: {pump['model']}**")
                st.caption(pump['spec'])
                
                with st.expander("📝 공학적 선정 조건 및 수식 근거 보기", expanded=False):
                    st.markdown(f"""
                    **[선정 근거 리포트]**
                    1. **유량 조건 ($Q$):** 시스템 전체 지배 노드의 총 유입량 {total_inflow:.3f} $m^3/s$에 대한 질량 보존 법칙 완벽 충족.
                    2. **손실 조건 ($\Delta P$):** 배관 선로의 지오메트리를 Hardy Cross 기법으로 수렴 연산한 결과 도출된 총 마찰 손실압 **{total_dp_loss:,.1f} $N/m^2$**을 안정적으로 밀어낼 수 있는 수두 능력 확보.
                    3. **동력 사양 ($W_{{pump}}$):** 이론 동력({required_power_kw:.2f} kW)에 기계적 효율($\\eta={pump_eff:.2f}$) 및 산업용 안전 마진(약 15%)을 가산하여 최적 수력 도달 기종으로 자동 매칭함.
                    
                    **[지배 방정식]**
                    $$W_{{pump}} = \\frac{{\\Delta P \\cdot Q}}{{\\eta \\cdot 1000}} \\quad [kW]$$
                    """)
                
                catalog_url = "https://product-selection.grundfos.com/?lc=KOR"
                st.link_button(
                    f"⚙️ {pump['search_name']} 계열 카탈로그에서 세부 사양서 매칭하기", 
                    catalog_url,
                    key=unique_key,
                    use_container_width=True
                )

    st.divider()
    
    st.subheader("🧐 열유체 공학적 설계 종합 진단 소견 (진동/소음 예측)")
    avg_diameter = sum(p.D for p in pipes) / len(pipes)
    
    col_eval1, col_eval2 = st.columns(2)
    with col_eval1:
        st.metric("배관망 평균 관경 (Avg D)", f"{avg_diameter:.3f} m")
    with col_eval2:
        if total_dp_loss < 50000:
            st.success(f"🎉 **설계 합격 (압력 최적화 달성):** 현재 전체 압력 손실치({total_dp_loss:,.1f} N/m²)가 경제적 안정 범위 내에 있습니다. 유체 마찰 저항이 효과적으로 억제되었습니다.")
        else:
            st.warning(f"⚠️ **압력 저하 설계 보완 필요:** 관로 손실 압력이 {total_dp_loss:,.1f} N/m²로 다소 높습니다. '전체 압력 저하'를 위해 직경(D)을 확장하는 설계 피드백을 권장합니다.")

    high_vibration_pipes = [p for p in pipes if p.get_delta_p_per_km() >= 557.0]
    if high_vibration_pipes:
        bad_ids = ", ".join([f"Pipe #{p.id}" for p in high_vibration_pipes])
        st.error(f"🚨 **배관 파괴 및 과도한 진동 경고:** 현재 계통 내 [{bad_ids}] 선로의 단위 압력 강하가 교재 소음 발생 기준치(557 kPa/km)를 초과했습니다! 장기 가동 시 허용 불가능한 공진 및 구조 결함이 예상되므로 관경 확장이 시급합니다.")
    else:
        st.info("✅ **진동/소음 안전성 검증:** 모든 배관 선로의 1000m당 압력 손실이 기준치(557 kPa/km) 미만으로 유지되어 매우 정숙하고 안정적인 유동 거동이 보장됩니다.")

    st.write(" ")
    st.markdown("### 📝 [최종 설계 출력 사양서 (Specification Summary)]")
    spec_data = {
        "설계 항목 (Design Item)": [
            "배관망 설계 구조 (Geometry Layout)", 
            "감지된 독립 루프개수 (Closed Loops)", 
            "계통 평균 관경 (Average Pipe Size)", 
            "총 마찰 손실 압력 (Total Head Loss)", 
            "공학적 권장 추천 펌프 라인업 (Recommended Series)",
            "최종 소요 전력 모터 출력 (Motor Power P2)"
        ],
        "최종 권장 사양 (Recommended Specifications)": [
            f"총 {len(pipes)}개 관로 분기망 구축",
            f"{len(auto_loops)}개 독립 폐회로 위상 제어",
            f"{avg_diameter*1000:.1f} mm (표준 스케일 실척 반영)",
            f"{total_dp_loss/1000:.2f} kPa",
            f"Grundfos {series_name.split(' ')[0]} Line업 (수력 범위 매칭)",
            f"최소 모터 정격 {(required_power_kw*1.15):.2f} kW 사양 권장 (안전율 15% 가산)"
        ]
    }
    st.table(pd.DataFrame(spec_data))

    st.subheader("📋 파이프 라인별 해석 결과 상세 내역")
    result_table = []
    for p in pipes:
        vib_status = "⚠️ 위험" if p.get_delta_p_per_km() >= 557.0 else "✅ 안전"
        result_table.append({
            "배관 번호": f"Pipe {p.id}",
            "연결 노드": f"{p.start} ➔ {p.end}",
            "최종 유량 (m³/s)": f"{p.Q:.4f}",
            "마찰 계수 (f)": f"{p.get_friction_factor():.4f}",
            "압력 강하 (N/m²)": f"{p.get_delta_p():,.1f}",
            "단위 압력 손실 (kPa/km)": f"{p.get_delta_p_per_km():.2f}",
            "진동/소음 예측": vib_status
        })
    st.dataframe(result_table, use_container_width=True)
