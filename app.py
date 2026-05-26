import streamlit as st
import economic_dia  # 분리한 경제지름 모듈 임포트

# --- 1. 페이지 공통 설정 (최상단에 한 번만 선언해야 합니다) ---
st.set_page_config(
    page_title="공학용 유체 설계 및 배관망 해석 시스템", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- 2. 사이드바 메인 제어 내비게이션 ---
st.sidebar.title("🛠️ 유체 시스템 플랫폼 v9.0")
st.sidebar.markdown("### Capstone Design Project")

# 프로젝트 선택 라디오 버튼
project_mode = st.sidebar.radio(
    "설계 및 해석 모드 선택",
    (
        "📈 Darcy-Weisbach 경제지름 구하기", 
        "🕸️ Hardy Cross 배관망 해석"
    )
)

st.sidebar.divider()
st.sidebar.caption("개발자: 이기원 (Mechanical Engineering)")
st.sidebar.caption("지도교수: 배주열 교수님")

# --- 3. 사용자가 선택한 모드에 따른 화면 분기 ---
if project_mode == "📈 Darcy-Weisbach 경제지름 구하기":
    # 메인 타이틀 출력
    st.title("🚀 공학용 유체 수송 설계 시스템")
    st.caption("기존 단일 배관 라인의 이론적 최적화 도출 및 상용 규격/비용 분석 통합 시스템")
    st.write("---")
    
    # economic_dia.py의 메인 함수 실행
    economic_dia.run_economic_dia()

elif project_mode == "🕸️ Hardy Cross 배관망 해석":
    # 신규 프로젝트 메인 타이틀 출력
    st.title("🕸️ Hardy Cross 배관망 해석 및 최적 펌프 선정 시스템")
    st.caption("다중 루프 배관망에서의 유량 배분, 압력 강하 시뮬레이션 및 민감도 기반 제어 노드 탐색")
    st.write("---")
    
    # 임시 안내 메시지 (이후 hardy_cross.py 연동 시 이 부분을 수정하게 됩니다)
    st.info("💡 현재 배관망 해석 엔진 및 NetworkX 기반 2D CAD 시각화 모듈을 준비 중입니다.")
    
    # 교수님 피드백 핵심 요약 가이드라인 배치 (개발 방향성 리마인드용)
    with st.expander("📌 교수님 피드백 및 개발 요구사항 (체크리스트)", expanded=True):
        st.markdown("""
        1. **기본 가정 사항**
            - 노드별 높이 차이($z$) 없음, 부차 손실(Minor Loss) 무시
            - 배관 지오메트리(길이, 직경) 고정 후 각 파이프의 유량 및 압력 강하($\\Delta p$) 수치 해석
        2. **Chen 방정식을 이용한 마찰 계수($f$) 업데이트**
            - 난류 영역에서의 정확한 수렴을 위해 매 Iteration마다 $Re$수에 따른 $f$값 자동 갱신
        3. **핵심 공학적 의사결정 기능 (민감도 분석)**
            - 특정 노드에서 유량을 컨트롤하는 것이 전체 시스템 밸런스에 가장 유리한지 판별
            - 펌프의 용량/유량을 어느 노드에서 변화시키는 것이 전체 압력 손실을 최소화하는지 탐색
        4. **최적 펌프 선정 연계**
            - 배관망의 최종 압력 강하 곡선(System Curve)을 도출하고, 이를 만족하는 정격 펌프 사양(동력, 양정) 자동 연산
        """)
        
    st.warning("👨‍💻 하디 크로스 백엔드 수식 로직과 2D 도면 UI 코딩 준비가 완료되면 이 창에 실시간 그래프와 함께 시뮬레이터가 구동됩니다.")
