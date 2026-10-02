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


# ============================================================
# ④ DNA 변이 분석
# ============================================================

st.divider()
st.header("④ DNA 변이 분석")

st.write(
    "기준 DNA와 비교 DNA의 염기서열을 비교하거나, "
    "ClinVar에서 알려진 유전자 변이를 검색할 수 있습니다."
)


# ------------------------------------------------------------
# NCBI Nuccore 검색
# ------------------------------------------------------------
def search_nuccore(query, retmax=10):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "nuccore",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "api_key": NCBI_API_KEY,
        "email": NCBI_EMAIL
    }

    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()

    data = response.json()

    ids = data["esearchresult"]["idlist"]
    count = int(data["esearchresult"]["count"])

    return ids, count


# ------------------------------------------------------------
# NCBI 서열 정보 가져오기
# ------------------------------------------------------------
def get_nuccore_summaries(ids):

    if not ids:
        return []

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    params = {
        "db": "nuccore",
        "id": ",".join(ids),
        "retmode": "json",
        "api_key": NCBI_API_KEY,
        "email": NCBI_EMAIL
    }

    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()

    data = response.json()

    results = []

    for uid in ids:

        info = data["result"].get(uid, {})

        results.append({
            "uid": uid,
            "title": info.get("title", "제목 없음"),
            "accession": info.get("accessionversion", ""),
            "length": info.get("slen", 0)
        })

    return results


# ------------------------------------------------------------
# NCBI에서 실제 DNA 서열 가져오기
# ------------------------------------------------------------
def fetch_nuccore_sequence(uid):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    params = {
        "db": "nuccore",
        "id": uid,
        "rettype": "fasta",
        "retmode": "text",
        "api_key": NCBI_API_KEY,
        "email": NCBI_EMAIL
    }

    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()

    fasta = response.text

    # FASTA 첫 줄(>로 시작하는 설명)을 제외하고
    # 실제 DNA 염기서열만 합침
    sequence = "".join(
        line.strip()
        for line in fasta.splitlines()
        if not line.startswith(">")
    )

    return sequence.upper()


# ------------------------------------------------------------
# 두 DNA 서열의 변이 탐색
# ------------------------------------------------------------
def find_variants(reference, sample):

    if len(reference) != len(sample):
        return None

    variants = []

    for i, (ref_base, sample_base) in enumerate(zip(reference, sample)):

        if ref_base != sample_base:

            variants.append({
                "위치": i + 1,
                "기준 염기": ref_base,
                "비교 염기": sample_base
            })

    return variants


# ------------------------------------------------------------
# Transition / Transversion 구분
# ------------------------------------------------------------
def classify_substitution(ref, alt):

    transition_pairs = {
        ("A", "G"),
        ("G", "A"),
        ("C", "T"),
        ("T", "C")
    }

    if (ref, alt) in transition_pairs:
        return "Transition"

    return "Transversion"


# ------------------------------------------------------------
# 변이가 코돈과 아미노산에 미치는 영향 분석
# ------------------------------------------------------------
def analyze_codon_change(reference, sample, position):

    # position은 1부터 시작하므로 Python index로 변환
    index = position - 1

    codon_start = (index // 3) * 3

    ref_codon = reference[codon_start:codon_start + 3]
    alt_codon = sample[codon_start:codon_start + 3]

    # 마지막에 불완전한 코돈이 있으면 분석하지 않음
    if len(ref_codon) != 3 or len(alt_codon) != 3:
        return None

    ref_aa = CODON_TABLE.get(ref_codon, "?")
    alt_aa = CODON_TABLE.get(alt_codon, "?")

    # 변이 종류 분류
    if ref_aa == alt_aa:
        effect = "동의적 변이 (Synonymous)"

    elif alt_aa == "*":
        effect = "종결 코돈 생성 (Nonsense)"

    else:
        effect = "아미노산 변화 (Missense)"

    return {
        "기준 코돈": ref_codon,
        "변이 코돈": alt_codon,
        "기준 아미노산": ref_aa,
        "변이 아미노산": alt_aa,
        "영향": effect
    }


# ============================================================
# 데이터 소스 선택
# ============================================================

source = st.radio(
    "비교 데이터 선택",
    [
        "직접 DNA 서열 입력",
        "NCBI에서 비교 서열 검색",
        "ClinVar에서 알려진 변이 검색"
    ],
    horizontal=True
)


# ============================================================
# 1. 직접 DNA 입력
# ============================================================

if source == "직접 DNA 서열 입력":

    st.subheader("비교 DNA 직접 입력")

    sample_input = st.text_area(
        "비교할 DNA 서열",
        height=150,
        placeholder="ATGCGT..."
    )

    if st.button(
        "DNA 비교 분석",
        type="primary",
        use_container_width=True,
        key="manual_variant"
    ):

        if "sequence" not in st.session_state:

            st.warning("먼저 NCBI에서 기준 DNA를 검색해주세요.")

        else:

            sample = (
                sample_input
                .upper()
                .replace(" ", "")
                .replace("\n", "")
            )

            if not sample:

                st.warning("비교할 DNA 서열을 입력해주세요.")

            elif any(base not in "ATGC" for base in sample):

                st.error("DNA 서열에는 A, T, G, C만 입력해주세요.")

            else:

                reference = st.session_state.sequence

                variants = find_variants(reference, sample)

                if variants is None:

                    st.error(
                        "현재 분석 방식에서는 기준 DNA와 비교 DNA의 "
                        "길이가 같아야 합니다."
                    )

                else:

                    st.session_state.variant_results = variants
                    st.session_state.variant_sample = sample


# ============================================================
# 2. NCBI에서 비교 DNA 검색
# ============================================================

elif source == "NCBI에서 비교 서열 검색":

    st.subheader("NCBI 비교 서열 검색")

    compare_query = st.text_input(
        "유전자명, accession 또는 검색어 입력",
        placeholder="예: TP53 Homo sapiens"
    )

    if st.button(
        "NCBI 비교 서열 검색",
        use_container_width=True
    ):

        if compare_query:

            try:

                ids, total_count = search_nuccore(compare_query)

                st.session_state.compare_ids = ids
                st.session_state.compare_total_count = total_count

                if ids:
                    st.session_state.compare_summaries = (
                        get_nuccore_summaries(ids)
                    )

            except Exception as e:

                st.error(f"NCBI 검색 중 오류가 발생했습니다: {e}")


    # 검색 결과 유지
    if "compare_summaries" in st.session_state:

        summaries = st.session_state.compare_summaries

        total = st.session_state.get(
            "compare_total_count",
            len(summaries)
        )

        st.success(
            f"NCBI에서 총 {total:,}개의 검색 결과를 찾았습니다. "
            f"상위 {len(summaries)}개를 표시합니다."
        )

        options = {}

        for item in summaries:

            label = (
                f"{item['accession']} | "
                f"{item['title']} | "
                f"{item['length']:,} bp"
            )

            options[label] = item["uid"]

        selected_label = st.selectbox(
            "비교할 서열 선택",
            list(options.keys())
        )

        if st.button(
            "선택한 서열과 비교",
            type="primary",
            use_container_width=True
        ):

            if "sequence" not in st.session_state:

                st.warning("먼저 기준 DNA를 검색해주세요.")

            else:

                try:

                    uid = options[selected_label]

                    sample = fetch_nuccore_sequence(uid)
                    reference = st.session_state.sequence

                    variants = find_variants(reference, sample)

                    if variants is None:

                        st.error(
                            "두 서열의 길이가 다릅니다. "
                            "현재 버전에서는 길이가 같은 서열끼리 "
                            "염기 위치를 직접 비교합니다."
                        )

                        st.info(
                            f"기준 서열: {len(reference):,} bp / "
                            f"비교 서열: {len(sample):,} bp"
                        )

                    else:

                        st.session_state.variant_results = variants
                        st.session_state.variant_sample = sample

                        st.success(
                            "NCBI 비교 서열을 불러왔습니다."
                        )

                except Exception as e:

                    st.error(
                        f"서열을 가져오는 중 오류가 발생했습니다: {e}"
                    )


# ============================================================
# 3. ClinVar 알려진 변이 검색
# ============================================================

elif source == "ClinVar에서 알려진 변이 검색":

    st.subheader("ClinVar 알려진 변이 검색")

    st.caption(
        "ClinVar에서는 DNA 전체 서열 대신 "
        "보고된 인간 유전체 변이 정보를 검색합니다."
    )

    clinvar_query = st.text_input(
        "유전자명 또는 변이 검색",
        placeholder="예: BRCA1 또는 TP53"
    )

    if st.button(
        "ClinVar 검색",
        type="primary",
        use_container_width=True
    ):

        if clinvar_query:

            try:

                # ClinVar 검색
                search_url = (
                    "https://eutils.ncbi.nlm.nih.gov/"
                    "entrez/eutils/esearch.fcgi"
                )

                search_params = {
                    "db": "clinvar",
                    "term": clinvar_query,
                    "retmode": "json",
                    "retmax": 10,
                    "api_key": NCBI_API_KEY,
                    "email": NCBI_EMAIL
                }

                response = requests.get(
                    search_url,
                    params=search_params,
                    timeout=20
                )

                response.raise_for_status()

                search_data = response.json()

                ids = search_data["esearchresult"]["idlist"]

                total_count = int(
                    search_data["esearchresult"]["count"]
                )

                st.session_state.clinvar_total = total_count


                if not ids:

                    st.warning("ClinVar 검색 결과가 없습니다.")

                    st.session_state.pop(
                        "clinvar_results",
                        None
                    )

                else:

                    # 검색된 ClinVar ID들의 상세 요약 가져오기
                    summary_url = (
                        "https://eutils.ncbi.nlm.nih.gov/"
                        "entrez/eutils/esummary.fcgi"
                    )

                    summary_params = {
                        "db": "clinvar",
                        "id": ",".join(ids),
                        "retmode": "json",
                        "api_key": NCBI_API_KEY,
                        "email": NCBI_EMAIL
                    }

                    summary_response = requests.get(
                        summary_url,
                        params=summary_params,
                        timeout=20
                    )

                    summary_response.raise_for_status()

                    summary_data = summary_response.json()

                    results = []

                    for uid in ids:

                        record = summary_data["result"].get(
                            uid,
                            {}
                        )

                        results.append({
                            "ClinVar ID": uid,
                            "변이": record.get(
                                "title",
                                "정보 없음"
                            ),
                            "변이 유형": record.get(
                                "obj_type",
                                record.get(
                                    "variation_set",
                                    "정보 없음"
                                )
                            ),
                            "임상적 분류": record.get(
                                "germline_classification",
                                {}
                            ).get(
                                "description",
                                "정보 없음"
                            )
                            if isinstance(
                                record.get(
                                    "germline_classification",
                                    {}
                                ),
                                dict
                            )
                            else "정보 없음"
                        })

                    st.session_state.clinvar_results = results

            except Exception as e:

                st.error(
                    f"ClinVar 검색 중 오류가 발생했습니다: {e}"
                )


    # ClinVar 결과를 rerun 후에도 유지
    if "clinvar_results" in st.session_state:

        results = st.session_state.clinvar_results

        total = st.session_state.get(
            "clinvar_total",
            len(results)
        )

        st.success(
            f"ClinVar에서 총 {total:,}개의 관련 변이를 찾았습니다. "
            f"상위 {len(results)}개를 표시합니다."
        )

        clinvar_df = pd.DataFrame(results)

        st.dataframe(
            clinvar_df,
            use_container_width=True,
            hide_index=True
        )

        st.info(
            "ClinVar 결과는 알려진 변이에 대한 데이터베이스 정보를 "
            "보여주는 기능입니다. 이 앱이 자체적으로 질병 여부를 "
            "판단하는 것은 아닙니다."
        )


# ============================================================
# DNA 직접 비교 결과
# ============================================================

if (
    source != "ClinVar에서 알려진 변이 검색"
    and "variant_results" in st.session_state
    and "variant_sample" in st.session_state
):

    variants = st.session_state.variant_results
    sample = st.session_state.variant_sample
    reference = st.session_state.sequence

    st.divider()
    st.subheader("변이 분석 결과")

    if len(variants) == 0:

        st.success(
            "두 DNA 서열 사이에서 염기 치환이 발견되지 않았습니다."
        )

    else:

        # ----------------------------------------------------
        # 기본 통계
        # ----------------------------------------------------

        transition_count = 0
        transversion_count = 0

        detailed_results = []

        for variant in variants:

            ref = variant["기준 염기"]
            alt = variant["비교 염기"]

            substitution_type = classify_substitution(
                ref,
                alt
            )

            if substitution_type == "Transition":
                transition_count += 1
            else:
                transversion_count += 1

            codon_info = analyze_codon_change(
                reference,
                sample,
                variant["위치"]
            )

            result = {
                "위치": variant["위치"],
                "기준 염기": ref,
                "비교 염기": alt,
                "치환 유형": substitution_type
            }

            if codon_info:

                result.update(codon_info)

            detailed_results.append(result)


        # ----------------------------------------------------
        # 결과 요약
        # ----------------------------------------------------

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                "총 변이 수",
                len(variants)
            )

        with col2:
            st.metric(
                "Transition",
                transition_count
            )

        with col3:
            st.metric(
                "Transversion",
                transversion_count
            )


        # ----------------------------------------------------
        # 상세 변이 표
        # ----------------------------------------------------

        st.subheader("염기 · 코돈 · 아미노산 변화")

        variant_df = pd.DataFrame(
            detailed_results
        )

        st.dataframe(
            variant_df,
            use_container_width=True,
            hide_index=True
        )


        # ----------------------------------------------------
        # 변이 위치 그래프
        # ----------------------------------------------------

        st.subheader("DNA 내 변이 위치")

        graph_df = pd.DataFrame({
            "위치": [
                v["위치"]
                for v in variants
            ],
            "변이": [
                f"{v['기준 염기']} → {v['비교 염기']}"
                for v in variants
            ],
            "값": [1] * len(variants)
        })

        fig = px.scatter(
            graph_df,
            x="위치",
            y="값",
            hover_name="변이",
            title="DNA 서열 내 변이 위치"
        )

        fig.update_yaxes(
            visible=False
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )
