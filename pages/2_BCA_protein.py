import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px


st.set_page_config(
    page_title="BCA 단백질 정량",
    page_icon="🧪"
)


# 표준용액 데이터를 이용해 회귀식을 계산하는 함수
def calculate_regression(concentrations, absorbances):

    x = np.array(concentrations)
    y = np.array(absorbances)

    slope, intercept = np.polyfit(x, y, 1)

    predicted = slope * x + intercept

    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    if ss_tot == 0:
        r_squared = 0
    else:
        r_squared = 1 - ss_res / ss_tot

    return slope, intercept, r_squared


# 회귀식으로 미지 시료 농도를 계산하는 함수
def calculate_unknown(absorbance, slope, intercept):

    if slope == 0:
        return None

    return (absorbance - intercept) / slope


# --------------------------------------------------
# 화면
# --------------------------------------------------

st.title("🧪 BCA 단백질 정량")

st.write(
    "BCA 실험에서 얻은 표준용액 데이터를 이용하여 "
    "미지 시료의 단백질 농도를 계산합니다."
)

st.subheader("표준용액 데이터")

default_data = pd.DataFrame({
    "단백질 농도 (mg/mL)": [
        0.0, 0.25, 0.5, 1.0, 1.5, 2.0
    ],
    "흡광도": [
        0.00, 0.10, 0.21, 0.41, 0.62, 0.82
    ]
})

standard_data = st.data_editor(
    default_data,
    num_rows="dynamic",
    use_container_width=True
)

st.subheader("미지 시료")

unknown_absorbance = st.number_input(
    "미지 시료 흡광도",
    min_value=0.0,
    value=0.30,
    step=0.01
)

if st.button(
    "BCA 분석하기",
    type="primary"
):

    data = standard_data.dropna()

    if len(data) < 2:

        st.error(
            "표준용액 데이터가 최소 2개 필요합니다."
        )

    else:

        concentrations = (
            data["단백질 농도 (mg/mL)"]
            .tolist()
        )

        absorbances = (
            data["흡광도"]
            .tolist()
        )

        if any(value < 0 for value in concentrations):

            st.error(
                "단백질 농도는 음수가 될 수 없습니다."
            )

        elif any(value < 0 for value in absorbances):

            st.error(
                "흡광도는 음수가 될 수 없습니다."
            )

        elif len(set(concentrations)) != len(concentrations):

            st.error(
                "동일한 단백질 농도가 중복되어 있습니다."
            )

        else:

            slope, intercept, r_squared = (
                calculate_regression(
                    concentrations,
                    absorbances
                )
            )

            unknown_concentration = (
                calculate_unknown(
                    unknown_absorbance,
                    slope,
                    intercept
                )
            )

            st.divider()

            st.subheader("분석 결과")

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric(
                    "미지 시료 흡광도",
                    f"{unknown_absorbance:.3f}"
                )

            with col2:
                st.metric(
                    "단백질 농도",
                    f"{unknown_concentration:.3f} mg/mL"
                )

            with col3:
                st.metric(
                    "R²",
                    f"{r_squared:.4f}"
                )

            st.subheader("BCA 표준곡선")

            x = np.array(concentrations)
            y = np.array(absorbances)

            graph_data = pd.DataFrame({
                "단백질 농도": x,
                "흡광도": y
            })

            fig = px.scatter(
                graph_data,
                x="단백질 농도",
                y="흡광도",
                title="BCA 단백질 표준곡선"
            )

            fig.add_scatter(
                x=x,
                y=slope * x + intercept,
                mode="lines",
                name="선형 회귀선"
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

            st.subheader("회귀식")

            st.code(
                f"y = {slope:.4f}x + {intercept:.4f}"
            )
