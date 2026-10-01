import streamlit as st
import pandas as pd
import plotly.express as px
import requests


# --------------------------------------------------
# 기본 설정
# --------------------------------------------------

st.set_page_config(
    page_title="DNA 분석",
    page_icon="🧬",
    layout="wide"
)


# --------------------------------------------------
# NCBI API 설정
# --------------------------------------------------

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Streamlit Secrets에서 API 정보를 가져옵니다.
NCBI_API_KEY = st.secrets["NCBI_API_KEY"]
NCBI_EMAIL = st.secrets["NCBI_EMAIL"]


# --------------------------------------------------
# DNA 분석 함수
# --------------------------------------------------

def count_bases(sequence):
    """DNA의 각 염기 개수를 계산합니다."""

    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }


def calculate_gc(sequence):
    """DNA의 GC 함량을 계산합니다."""

    if len(sequence) == 0:
        return 0

    gc_count = sequence.count("G") + sequence.count("C")

    return gc_count / len(sequence) * 100


# --------------------------------------------------
# NCBI 검색
# --------------------------------------------------

def search_ncbi(gene_name, organism):
    """NCBI에서 입력한 유전자를 검색합니다."""

    search_term = (
        f"{gene_name}[Gene Name] "
        f"AND {organism}[Organism]"
    )

    params = {
        "db": "nuccore",
        "term": search_term,
        "retmode": "json",
        "retmax": 5,
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_URL + "esearch.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    data = response.json()

    return data["esearchresult"]["idlist"]


def get_ncbi_summary(ids):
    """검색된 DNA 기록의 정보를 가져옵니다."""

    params = {
        "db": "nuccore",
        "id": ",".join(ids),
        "retmode": "json",
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_URL + "esummary.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    data = response.json()["result"]

    results = []

    for record_id in ids:

        record = data.get(record_id, {})

        results.append({
            "ID": record_id,
            "제목": record.get(
                "title",
                "정보 없음"
            ),
            "길이": record.get(
                "slen",
                "정보 없음"
            )
        })

    return results


def get_ncbi_sequence(record_id):
    """선택한 NCBI 기록에서 DNA 서열을 가져옵니다."""

    params = {
        "db": "nuccore",
        "id": record_id,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "BioLabAnalyzer",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }

    response = requests.get(
        NCBI_URL + "efetch.fcgi",
        params=params,
        timeout=15
    )

    response.raise_for_status()

    lines = response.text.strip().splitlines()

    # FASTA의 설명 부분을 제외하고 DNA만 가져옵니다.
    sequence = "".join(
        line.strip()
        for line in lines
        if not line.startswith(">")
    )

    return sequence


# --------------------------------------------------
# 코돈 → 아미노산 변환표
# --------------------------------------------------

codon_table = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",

    "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S",

    "TAT": "Y", "TAC": "Y",

    "TAA": "*", "TAG": "*",

    "TGT": "C", "TGC": "C",
    "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L",

    "CCT": "P", "CCC": "P",
    "CCA": "P", "CCG": "P",

    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",

    "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I",
    "ATA": "I", "ATG": "M",

    "ACT": "T", "ACC": "T",
    "ACA": "T", "ACG": "T",

    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",

    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V",

    "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A",

    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",

    "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G"
}


amino_acid_names = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine"
}


def translate_dna(sequence):
    """DNA를 3개씩 읽어 아미노산 서열로 변환합니다."""

    protein = []

    for i in range(0, len(sequence) - 2, 3):

        codon = sequence[i:i + 3]

        protein.append(
            codon_table.get(codon, "?")
        )

    return "".join(protein)


def find_orfs(sequence):
    """3개의 reading frame에서 ORF를 찾습니다."""

    orfs = []

    stop_codons = {
        "TAA",
        "TAG",
        "TGA"
    }

    # DNA를 0, 1, 2번째 위치부터 각각 읽어봅니다.
    for frame in range(3):

        i = frame

        while i <= len(sequence) - 3:

            codon = sequence[i:i + 3]

            # 시작 코돈을 찾습니다.
            if codon == "ATG":

                start = i
                j = i + 3

                # 시작 코돈 이후의 종결 코돈을 찾습니다.
                while j <= len(sequence) - 3:

                    current = sequence[j:j + 3]

                    if current in stop_codons:

                        end = j + 3

                        dna = sequence[start:end]

                        protein = translate_dna(dna)

                        orfs.append({
                            "frame": frame + 1,
                            "start": start + 1,
                            "end": end,
                            "dna": dna,
                            "protein": protein
                        })

                        break

                    j += 3

            i += 3

    return orfs


# --------------------------------------------------
# 화면
# --------------------------------------------------

st.title("🧬 DNA Analysis")

st.write(
    "NCBI에서 실제 유전자 서열을 가져와 "
    "DNA의 특성과 ORF를 분석합니다."
)


# ==================================================
# ① NCBI 유전자 검색
# ==================================================

st.subheader("① NCBI 유전자 검색")

col1, col2 = st.columns(2)

with col1:

    gene_name = st.text_input(
        "유전자명",
        placeholder="예: TP53"
    )

with col2:

    organism = st.text_input(
        "생물 종",
        value="Homo sapiens"
    )


if st.button(
    "NCBI에서 검색하기",
    type="primary",
    use_container_width=True
):

    if not gene_name:

        st.warning(
            "검색할 유전자명을 입력해주세요."
        )

    else:

        try:

            with st.spinner(
                "NCBI에서 검색 중..."
            ):

                ids = search_ncbi(
                    gene_name,
                    organism
                )

                if not ids:

                    st.warning(
                        "검색 결과가 없습니다."
                    )

                else:

                    results = get_ncbi_summary(ids)

                    st.session_state.ncbi_results = results

        except requests.exceptions.RequestException:

            st.error(
                "NCBI 서버와 연결할 수 없습니다."
            )


# 검색 결과 표시
if "ncbi_results" in st.session_state:

    st.write("검색 결과")

    result_df = pd.DataFrame(
        st.session_state.ncbi_results
    )

    st.dataframe(
        result_df,
        use_container_width=True,
        hide_index=True
    )

    selected = st.selectbox(
        "분석할 서열을 선택하세요.",
        st.session_state.ncbi_results,
        format_func=lambda x:
            f"{x['ID']} | {x['제목']}"
    )

    if st.button(
        "이 서열 분석하기",
        use_container_width=True
    ):

        try:

            with st.spinner(
                "DNA 서열을 가져오는 중..."
            ):

                sequence = get_ncbi_sequence(
                    selected["ID"]
                )

                st.session_state.sequence = sequence

                # 새로운 DNA를 가져오면 이전 ORF 분석 결과를 초기화합니다.
                st.session_state.pop("orfs", None)

                st.success(
                    "DNA 서열을 가져왔습니다."
                )

        except requests.exceptions.RequestException:

            st.error(
                "DNA 서열을 가져오는 데 실패했습니다."
            )


# ==================================================
# ② DNA 기본 분석
# ==================================================

if "sequence" in st.session_state:

    sequence = st.session_state.sequence

    st.divider()

    st.subheader("② DNA 기본 분석")

    col1, col2 = st.columns(2)

    with col1:

        st.metric(
            "DNA 길이",
            f"{len(sequence):,} bp"
        )

    with col2:

        gc_content = calculate_gc(
            sequence
        )

        st.metric(
            "GC 함량",
            f"{gc_content:.2f}%"
        )

    counts = count_bases(
        sequence
    )

    base_df = pd.DataFrame({
        "염기": list(counts.keys()),
        "개수": list(counts.values())
    })

    fig = px.bar(
        base_df,
        x="염기",
        y="개수",
        title="DNA 염기 조성"
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )


# ==================================================
# ③ ORF 및 아미노산 분석
# ==================================================

st.divider()

st.subheader(
    "③ ORF 및 아미노산 분석"
)

st.write(
    "DNA에서 시작 코돈 ATG부터 "
    "종결 코돈까지 이어지는 ORF를 찾고, "
    "이를 아미노산 서열로 변환합니다."
)


# --------------------------------------------------
# ORF 분석 버튼
# --------------------------------------------------

if st.button(
    "ORF 분석하기",
    type="primary",
    use_container_width=True
):

    # ORF 분석 결과를 저장합니다.
    orfs = find_orfs(sequence)

    st.session_state.orfs = orfs


# --------------------------------------------------
# ORF 분석 결과
# --------------------------------------------------

# 버튼을 다시 누르지 않아도 저장된 결과를 보여줍니다.
if "orfs" in st.session_state:

    orfs = st.session_state.orfs

    if not orfs:

        st.warning(
            "완전한 ORF를 찾지 못했습니다."
        )

    else:

        st.success(
            f"{len(orfs)}개의 ORF를 찾았습니다."
        )


        # ------------------------------------------
        # ORF 선택
        # ------------------------------------------

        orf_options = {
            f"ORF {i + 1} | "
            f"Frame {orf['frame']} | "
            f"{orf['start']}–{orf['end']}": i
            for i, orf in enumerate(orfs)
        }

        selected_orf_label = st.selectbox(
            "분석할 ORF를 선택하세요.",
            list(orf_options.keys())
        )

        selected_index = orf_options[
            selected_orf_label
        ]

        selected_orf = orfs[
            selected_index
        ]


        # ------------------------------------------
        # 선택한 ORF 기본 정보
        # ------------------------------------------

        col1, col2 = st.columns(2)

        with col1:

            st.metric(
                "ORF 길이",
                f"{len(selected_orf['dna'])} bp"
            )

        with col2:

            st.metric(
                "아미노산 길이",
                f"{len(selected_orf['protein']) - 1} aa"
            )


        # ------------------------------------------
        # 아미노산 서열
        # ------------------------------------------

        st.write(
            "번역된 아미노산 서열"
        )

        st.code(
            selected_orf["protein"]
        )


        # ------------------------------------------
        # 특정 아미노산 검색
        # ------------------------------------------

        st.write(
            "아미노산 서열에서 원하는 아미노산을 찾습니다."
        )

        amino_options = {
            f"{name} ({code})": code
            for code, name
            in amino_acid_names.items()
        }

        selected_amino_name = st.selectbox(
            "검색할 아미노산",
            list(amino_options.keys())
        )

        selected_amino = amino_options[
            selected_amino_name
        ]


        # ------------------------------------------
        # 아미노산 위치 검색
        # ------------------------------------------

        protein = selected_orf[
            "protein"
        ]

        positions = [
            i + 1
            for i, amino in enumerate(protein)
            if amino == selected_amino
        ]


        if positions:

            st.success(
                f"{selected_amino_name} "
                f"{len(positions)}개 발견"
            )

            st.write(
                "아미노산 위치:",
                ", ".join(
                    map(str, positions)
                )
            )

        else:

            st.info(
                "해당 아미노산을 찾지 못했습니다."
            )
