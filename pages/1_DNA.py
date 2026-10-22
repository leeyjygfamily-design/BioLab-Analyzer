import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import requests
from Bio.Align import PairwiseAligner

st.set_page_config(page_title="DNA 분석", page_icon="🧬", layout="wide")

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API = {"tool": "BioLabAnalyzer", "email": st.secrets["NCBI_EMAIL"], "api_key": st.secrets["NCBI_API_KEY"]}

# -------------------- 공통 함수 --------------------
def ncbi_get(endpoint, db, **params):
    r = requests.get(NCBI_URL + endpoint, params={"db": db, **params, **API}, timeout=20)
    r.raise_for_status()
    return r


def search_db(db, term, retmax=10):
    data = ncbi_get("esearch.fcgi", db, term=term, retmode="json", retmax=retmax).json()["esearchresult"]
    return data["idlist"], int(data["count"])


def summaries(db, ids):
    if not ids:
        return []
    data = ncbi_get("esummary.fcgi", db, id=",".join(ids), retmode="json").json()["result"]
    return [data.get(i, {}) for i in ids]


def fetch_sequence(record_id):
    text = ncbi_get("efetch.fcgi", "nuccore", id=record_id, rettype="fasta", retmode="text").text
    return "".join(x.strip() for x in text.splitlines() if not x.startswith(">" )).upper()


def clean_dna(seq):
    return "".join(seq.upper().split())


def base_counts(seq):
    return {b: seq.count(b) for b in "ATGC"}


def gc_percent(seq):
    return 100 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0


# -------------------- 번역 / ORF --------------------
codon_table = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

amino_names = {
    "A":"Alanine","R":"Arginine","N":"Asparagine","D":"Aspartic acid","C":"Cysteine",
    "E":"Glutamic acid","Q":"Glutamine","G":"Glycine","H":"Histidine","I":"Isoleucine",
    "L":"Leucine","K":"Lysine","M":"Methionine","F":"Phenylalanine","P":"Proline",
    "S":"Serine","T":"Threonine","W":"Tryptophan","Y":"Tyrosine","V":"Valine"
}


def translate(seq):
    return "".join(codon_table.get(seq[i:i+3], "?") for i in range(0, len(seq)-2, 3))


def find_orfs(seq):
    orfs, stops = [], {"TAA", "TAG", "TGA"}
    for frame in range(3):
        for start in range(frame, len(seq)-2, 3):
            if seq[start:start+3] != "ATG":
                continue
            for end in range(start+3, len(seq)-2, 3):
                if seq[end:end+3] in stops:
                    dna = seq[start:end+3]
                    orfs.append({"frame": frame+1, "start": start+1, "end": end+3,
                                 "dna": dna, "protein": translate(dna)})
                    break
    return orfs


# -------------------- 서열 정렬 / 변이 영향 분석 --------------------
def align_variants(reference, sample):
    """기준 DNA에서 비교 DNA와 유사한 구간을 찾아 정렬합니다(정방향)."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    alignment = aligner.align(reference, sample)[0]
    coords = alignment.coordinates

    variants, mapping, insertions = [], {}, {}
    matches = columns = 0
    for k in range(coords.shape[1] - 1):
        r1, r2 = map(int, coords[0, k:k+2])
        s1, s2 = map(int, coords[1, k:k+2])
        if r2 > r1 and s2 > s1:
            for ri, si in zip(range(r1, r2), range(s1, s2)):
                ref, alt = reference[ri], sample[si]
                mapping[ri] = alt
                columns += 1
                if ref == alt:
                    matches += 1
                else:
                    variants.append({"기준 위치": ri+1, "기준 염기": ref,
                                     "비교 염기": alt, "유형": "SNP"})
        elif r2 > r1:
            for ri in range(r1, r2):
                mapping[ri] = ""
                columns += 1
                variants.append({"기준 위치": ri+1, "기준 염기": reference[ri],
                                 "비교 염기": "-", "유형": "Deletion"})
        elif s2 > s1:
            insertion = sample[s1:s2]
            insertions[r1] = insertions.get(r1, "") + insertion
            columns += len(insertion)
            variants.append({"기준 위치": r1, "기준 염기": "-",
                             "비교 염기": insertion, "유형": "Insertion"})

    return {
        "variants": variants, "mapping": mapping, "insertions": insertions,
        "start": int(coords[0, 0])+1, "end": int(coords[0, -1]),
        "identity": 100*matches/columns if columns else 0,
        "coverage": 100*(int(coords[1, -1])-int(coords[1, 0]))/len(sample),
        "aligned_bases": int(coords[1, -1])-int(coords[1, 0]),
    }


def substitution_type(ref, alt):
    if ref == "-" or alt == "-":
        return "-"
    return "Transition" if (ref, alt) in {
        ("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")
    } else "Transversion"


def analyze_orf_effect(reference, result, orf):
    """선택한 ORF 전체가 정렬되었을 때만 단백질 수준의 영향을 계산합니다."""
    if not orf:
        return None, "ORF 분석에서 ORF를 선택하면 단백질 영향도 확인할 수 있습니다."

    start, end = orf["start"]-1, orf["end"]  # Python 인덱스
    mapping, insertions = result["mapping"], result["insertions"]
    if any(i not in mapping for i in range(start, end)):
        return None, "선택한 ORF 전체가 비교 서열과 정렬되지 않아 단백질 영향을 확정할 수 없습니다."

    ref_dna = reference[start:end]
    pieces = []
    for i in range(start, end):
        # 삽입 위치는 기준 염기의 바로 앞(0-based 경계)
        pieces.append(insertions.get(i, "") if i > start else "")
        pieces.append(mapping[i])
    alt_dna = "".join(pieces)

    relevant = [v for v in result["variants"]
                if (start+1 <= v["기준 위치"] <= end if v["유형"] != "Insertion"
                    else start < v["기준 위치"] < end)]
    if not relevant:
        return {"rows": [], "ref_protein": translate(ref_dna),
                "alt_protein": translate(alt_dna), "summary": "선택한 ORF에서 변이가 발견되지 않았습니다."}, None

    indels = [v for v in relevant if v["유형"] != "SNP"]
    net_change = len(alt_dna)-len(ref_dna)
    if indels:
        label = ("Frameshift 가능성 (삽입·결실의 순길이가 3의 배수가 아님)"
                 if net_change % 3 else
                 "In-frame 삽입·결실 가능성 (순길이가 3의 배수)")
        return {"rows": [], "ref_protein": translate(ref_dna),
                "alt_protein": translate(alt_dna), "summary": label,
                "indels": True}, None

    # 한 코돈에 SNP가 여러 개 있어도 최종 변이 코돈으로 함께 평가
    rows = []
    for codon_index in sorted({(v["기준 위치"]-1-start)//3 for v in relevant}):
        j = codon_index*3
        before, after = ref_dna[j:j+3], alt_dna[j:j+3]
        aa1, aa2 = translate(before), translate(after)
        if aa1 == aa2:
            effect = "Synonymous (아미노산 유지)"
        elif aa2 == "*" and aa1 != "*":
            effect = "Nonsense (조기 종결 코돈)"
        elif aa1 == "*" and aa2 != "*":
            effect = "Stop-loss (종결 코돈 소실)"
        else:
            effect = "Missense (아미노산 치환)"
        rows.append({"코돈 위치(ORF 기준)": codon_index+1,
                     "기준 코돈": before, "변이 코돈": after,
                     "기준 아미노산": aa1, "변이 아미노산": aa2,
                     "예측 영향": effect})
    return {"rows": rows, "ref_protein": translate(ref_dna),
            "alt_protein": translate(alt_dna),
            "summary": "선택한 ORF의 코돈·아미노산 변화를 분석했습니다."}, None


def recommend_orfs(orfs, result):
    """정렬된 기준 좌표로 ORF 전체 포함 여부와 겹침 비율을 평가합니다."""
    covered = result["mapping"]
    ranked = []
    for idx, orf in enumerate(orfs):
        start, end = orf["start"] - 1, orf["end"]
        length = end - start
        overlap = sum(i in covered for i in range(start, end))
        if not overlap:
            continue
        fraction = overlap / length
        ranked.append({"index": idx, "orf": orf, "fraction": fraction,
                       "full": overlap == length, "overlap": overlap,
                       "length_diff": abs(length - len(st.session_state.get("variant_sample", "")))})
    ranked.sort(key=lambda x: (-int(x["full"]), -x["fraction"],
                               -x["overlap"], x["length_diff"]))
    return ranked


def dna_variant_figure(start, end, variants):
    """정렬된 DNA 구간을 이중나선 모식도로 표현합니다. 실제 3차원 구조는 아닙니다."""
    span = max(1, end - start)
    x = np.linspace(start, end, 500)
    theta = 4 * np.pi * (x - start) / span
    fig = go.Figure()
    for phase, color, label in [(0, "#2979B8", "DNA 가닥 1"),
                                (np.pi, "#D18A35", "DNA 가닥 2")]:
        fig.add_trace(go.Scatter(x=x, y=np.sin(theta + phase), mode="lines",
                                 line=dict(color=color, width=3), name=label,
                                 hoverinfo="skip"))
    for xi in np.linspace(start, end, 32):
        yy = np.sin(4 * np.pi * (xi - start) / span)
        fig.add_shape(type="line", x0=xi, x1=xi, y0=-yy, y1=yy,
                      line=dict(color="rgba(130,130,130,0.35)", width=1))
    colors = {"SNP": "#DC3545", "Deletion": "#7B4DC4", "Insertion": "#159A75"}
    for typ in ("SNP", "Deletion", "Insertion"):
        group = [v for v in variants if v["유형"] == typ]
        if not group:
            continue
        fig.add_trace(go.Scatter(
            x=[v["기준 위치"] for v in group], y=[0]*len(group),
            mode="markers", name=typ,
            marker=dict(size=12, color=colors[typ], symbol="diamond",
                        line=dict(color="white", width=1)),
            customdata=[[v["기준 염기"], v["비교 염기"]] for v in group],
            hovertemplate="위치 %{x}<br>기준 %{customdata[0]} → 비교 %{customdata[1]}<extra>"+typ+"</extra>"
        ))
    fig.update_layout(title="DNA 이중나선 모식도와 변이 위치", xaxis_title="기준 DNA 위치 (bp)",
                      yaxis=dict(visible=False, range=[-1.35, 1.35]),
                      height=350, legend=dict(orientation="h", y=-0.22),
                      margin=dict(l=10, r=10, t=55, b=70))
    return fig


def save_alignment(sample):
    reference = st.session_state.sequence
    # 로컬 정렬은 큰 서열에서 느리거나 메모리를 많이 사용할 수 있음
    if len(reference) > 30000 or len(sample) > 10000:
        st.warning("현재 버전은 기준 30,000 bp, 비교 10,000 bp 이하만 지원합니다.")
        return
    try:
        result = align_variants(reference, sample)
        if result["aligned_bases"] < 12 or result["coverage"] < 50:
            st.warning("유사한 구간이 너무 짧거나 비교 서열의 절반 미만만 정렬되었습니다. 결과를 저장하지 않았습니다.")
            return
        st.session_state.pop("recommended_orf_choice", None)
        st.session_state.variant_result = result
        st.session_state.variant_sample = sample
    except (ValueError, IndexError, MemoryError) as e:
        st.error(f"서열 정렬에 실패했습니다: {e}")


# ==================== 화면 ====================
st.title("🧬 DNA Analysis")
st.write("NCBI 유전자 서열을 가져와 DNA 특성, ORF, 아미노산 및 변이를 분석합니다.")

# ① NCBI 검색
st.subheader("① NCBI 유전자 검색")
c1, c2 = st.columns(2)
gene = c1.text_input("유전자명", placeholder="예: TP53")
organism = c2.text_input("생물 종", value="Homo sapiens")

if st.button("NCBI에서 검색하기", type="primary", use_container_width=True):
    if not gene:
        st.warning("검색할 유전자명을 입력해주세요.")
    else:
        try:
            ids, total = search_db("nuccore", f"{gene}[Gene Name] AND {organism}[Organism]", 5)
            rows = summaries("nuccore", ids)
            st.session_state.ncbi_total = total
            st.session_state.ncbi_results = [
                {"ID": i, "제목": r.get("title", "정보 없음"), "길이": r.get("slen", "정보 없음")}
                for i, r in zip(ids, rows)
            ]
        except requests.RequestException:
            st.error("NCBI 서버와 연결할 수 없습니다.")

if "ncbi_results" in st.session_state:
    results = st.session_state.ncbi_results
    st.write(f"검색 결과: 총 {st.session_state.get('ncbi_total', len(results)):,}개 (상위 {len(results)}개 표시)")
    st.dataframe(pd.DataFrame(results), use_container_width=True, hide_index=True)
    selected = st.selectbox("분석할 서열", results, format_func=lambda x: f"{x['ID']} | {x['제목']}")

    if st.button("이 서열 분석하기", use_container_width=True):
        try:
            st.session_state.sequence = fetch_sequence(selected["ID"])
            for key in ("orfs", "selected_orf", "variant_result", "variant_sample", "recommended_orf_choice"):
                st.session_state.pop(key, None)
            st.success("DNA 서열을 가져왔습니다.")
        except requests.RequestException:
            st.error("DNA 서열을 가져오는 데 실패했습니다.")


if "sequence" not in st.session_state:
    st.info("먼저 위에서 기준 DNA 서열을 선택해주세요.")
    st.stop()

sequence = st.session_state.sequence

# ② 기본 분석
st.divider()
st.subheader("② DNA 기본 분석")
c1, c2 = st.columns(2)
c1.metric("DNA 길이", f"{len(sequence):,} bp")
c2.metric("GC 함량", f"{gc_percent(sequence):.2f}%")
counts = base_counts(sequence)
fig = px.bar(pd.DataFrame({"염기": counts.keys(), "개수": counts.values()}), x="염기", y="개수", title="DNA 염기 조성")
st.plotly_chart(fig, use_container_width=True)
with st.expander("GC 함량은 무엇을 의미하나요?"):
    st.markdown("""**GC 함량(%) = (G의 개수 + C의 개수) ÷ 전체 염기 수 × 100**

DNA에서 구아닌(G)과 사이토신(C)이 차지하는 비율을 의미합니다. G-C 염기쌍은 수소 결합 3개, A-T 염기쌍은 2개를 형성하므로, **다른 조건이 같다면 GC 함량이 높은 이중가닥 DNA는 일반적으로 녹는점(Tm)이 높아지는 경향**이 있습니다.

GC 함량은 PCR 프라이머 설계, 서열 특성 비교, DNA의 열적 안정성 추정에 활용됩니다. 다만 **GC 함량만으로 DNA의 실제 안정성이나 유전자 발현량, 단백질 기능을 판단할 수는 없습니다.** DNA 길이·염기 배열·염 농도 등도 영향을 줍니다.""")

# ③ ORF
st.divider()
st.subheader("③ ORF 및 아미노산 분석")
with st.expander("ORF란?"):
    st.write("ORF(Open Reading Frame)는 시작 코돈에서 종결 코돈까지 이어져 단백질로 번역될 가능성이 있는 염기서열 구간입니다.")

if st.button("ORF 분석하기", type="primary", use_container_width=True):
    st.session_state.orfs = find_orfs(sequence)
    st.session_state.pop("selected_orf", None)
    st.session_state.pop("recommended_orf_choice", None)

selected_orf = None
if "orfs" in st.session_state:
    orfs = st.session_state.orfs
    if not orfs:
        st.warning("완전한 ORF를 찾지 못했습니다.")
    else:
        st.success(f"총 {len(orfs)}개의 ORF를 찾았습니다.")
        csort, cmin = st.columns(2)
        sort_method = csort.selectbox(
            "ORF 정렬 기준",
            ["긴 길이순", "짧은 길이순", "DNA 위치순", "비교 서열 길이에 가까운 순"]
        )
        min_length = cmin.number_input("최소 ORF 길이 (bp)", min_value=0, value=300, step=30)
        target_length = None
        if sort_method == "비교 서열 길이에 가까운 순":
            default_target = len(st.session_state.get("variant_sample", ""))
            target_length = st.number_input(
                "비교 DNA 길이 (bp)", min_value=1,
                value=max(1, default_target), step=1,
                help="비교 DNA 길이를 입력하면 길이 차이가 작은 ORF부터 표시합니다."
            )
        candidates = [(i, o) for i, o in enumerate(orfs) if len(o["dna"]) >= min_length]
        if sort_method == "긴 길이순":
            candidates.sort(key=lambda item: len(item[1]["dna"]), reverse=True)
        elif sort_method == "짧은 길이순":
            candidates.sort(key=lambda item: len(item[1]["dna"]))
        elif sort_method == "비교 서열 길이에 가까운 순":
            candidates.sort(key=lambda item: abs(len(item[1]["dna"]) - target_length))
        else:
            candidates.sort(key=lambda item: item[1]["start"])

        st.caption(f"조건에 맞는 ORF: {len(candidates)}개 / 전체 {len(orfs)}개")
        if not candidates:
            st.warning("조건에 맞는 ORF가 없습니다. 최소 길이를 낮춰보세요.")
            st.session_state.pop("selected_orf", None)
            st.stop()

        selected_id = st.selectbox(
            "분석할 ORF",
            [i for i, _ in candidates],
            format_func=lambda i: (
                f"ORF {i+1} | {len(orfs[i]['dna']):,} bp | "
                f"Frame {orfs[i]['frame']} | {orfs[i]['start']:,}–{orfs[i]['end']:,}"
            )
        )
        selected_orf = orfs[selected_id]
        st.session_state.selected_orf = selected_orf

        c1, c2 = st.columns(2)
        c1.metric("ORF 길이", f"{len(selected_orf['dna'])} bp")
        c2.metric("아미노산 길이", f"{max(0, len(selected_orf['protein'])-1)} aa")
        st.write("번역된 아미노산 서열")
        st.code(selected_orf["protein"])

        options = {f"{name} ({code})": code for code, name in amino_names.items()}
        amino_label = st.selectbox("검색할 아미노산", options)
        code = options[amino_label]
        positions = [i+1 for i, aa in enumerate(selected_orf["protein"]) if aa == code]
        st.write(f"{amino_label}: {len(positions)}개")
        if positions:
            st.write("위치:", ", ".join(map(str, positions)))
else:
    selected_orf = st.session_state.get("selected_orf")

# ④ 변이 분석
st.divider()
st.subheader("④ DNA 변이 및 단백질 영향 분석")
st.write("기준 DNA에서 비교 서열과 가장 유사한 영역을 정렬하고, 선택한 ORF의 코돈·아미노산 변화를 살펴봅니다.")
st.caption("정방향 국소 정렬을 사용합니다. 유사성이 낮거나 역방향인 서열은 별도 확인이 필요합니다.")

source = st.radio("비교 데이터", ["직접 DNA 입력", "NCBI 비교 서열"], horizontal=True)

if source == "직접 DNA 입력":
    sample_text = st.text_area("비교할 DNA 서열", height=130, placeholder="ATGCGT...")
    if st.button("변이 분석하기", type="primary", use_container_width=True):
        sample = clean_dna(sample_text)
        if not sample:
            st.warning("비교할 DNA를 입력해주세요.")
        elif any(b not in "ATGC" for b in sample):
            st.error("A, T, G, C만 입력해주세요.")
        else:
            save_alignment(sample)

else:
    query = st.text_input("유전자명, accession 또는 검색어", placeholder="예: TP53 Homo sapiens")
    if st.button("NCBI 비교 서열 검색", use_container_width=True):
        if not query.strip():
            st.warning("검색어를 입력해주세요.")
        else:
            try:
                ids, total = search_db("nuccore", query, 10)
                rows = summaries("nuccore", ids)
                st.session_state.compare_total = total
                st.session_state.compare_results = [
                    {"ID": i, "제목": r.get("title", "제목 없음"),
                     "Accession": r.get("accessionversion", ""),
                     "길이": r.get("slen", 0)}
                    for i, r in zip(ids, rows)
                ]
            except (requests.RequestException, KeyError, ValueError):
                st.error("NCBI 검색에 실패했습니다.")

    if "compare_results" in st.session_state:
        rows = st.session_state.compare_results
        st.write(f"검색 결과: 총 {st.session_state.get('compare_total', len(rows)):,}개 (상위 {len(rows)}개 표시)")
        if rows:
            choice = st.selectbox(
                "비교할 서열", rows,
                format_func=lambda x: f"{x['Accession']} | {x['제목']} | {x['길이']} bp"
            )
            if st.button("선택한 서열 변이 분석하기", type="primary", use_container_width=True):
                try:
                    sample = fetch_sequence(choice["ID"])
                    if not sample or any(b not in "ATGC" for b in sample):
                        st.error("선택한 기록에 A/T/G/C 이외의 문자가 있어 분석할 수 없습니다.")
                    else:
                        save_alignment(sample)
                except requests.RequestException:
                    st.error("비교 서열을 가져오지 못했습니다.")
        else:
            st.info("검색 결과가 없습니다.")

if "variant_result" in st.session_state:
    r = st.session_state.variant_result
    variants = r["variants"]
    st.divider()
    st.subheader("서열 정렬 결과")
    c1, c2, c3 = st.columns(3)
    c1.metric("기준 DNA의 유사 구간", f"{r['start']:,}–{r['end']:,} bp")
    c2.metric("정렬 구간 일치도", f"{r['identity']:.2f}%")
    c3.metric("비교 서열 정렬 비율", f"{r['coverage']:.2f}%")
    st.caption("일치도는 정렬된 열 기준이며, 정렬되지 않은 비교 서열의 끝부분은 변이로 집계하지 않습니다.")

    if variants:
        df = pd.DataFrame([
            {**v, "치환 유형": substitution_type(v["기준 염기"], v["비교 염기"])}
            for v in variants
        ])
        st.write(f"정렬된 영역에서 총 {len(df)}개 변이 항목을 찾았습니다.")
        st.dataframe(df, use_container_width=True, hide_index=True)
        counts = df["유형"].value_counts().rename_axis("유형").reset_index(name="개수")
        st.plotly_chart(px.bar(counts, x="유형", y="개수", title="변이 유형"), use_container_width=True)
        with st.expander("변이 유형 설명: SNP / Deletion / Insertion"):
            st.markdown("""- **SNP(단일 염기 치환):** 한 위치의 염기가 다른 염기로 바뀐 경우입니다. 코돈에 따라 아미노산이 유지되거나 달라질 수 있습니다.
- **Deletion(결실):** 기준 서열에 있던 염기가 비교 서열에서 사라진 경우입니다.
- **Insertion(삽입):** 기준 서열에는 없던 염기가 비교 서열에 추가된 경우입니다.

코딩 영역에서 3의 배수가 아닌 염기가 삽입·결실되면 읽기 틀이 달라질 수 있습니다.""")
        with st.expander("치환 유형 설명: Transition / Transversion"):
            st.markdown("""- **Transition(전이):** 같은 계열 염기끼리의 치환입니다. 퓨린 A↔G 또는 피리미딘 C↔T.
- **Transversion(전환):** 서로 다른 계열 간 치환입니다. A/G ↔ C/T.

이 분류는 **염기의 화학적 계열**을 나타내며, 그 자체로 단백질 기능에 미치는 영향의 크기를 뜻하지는 않습니다.""")
        st.plotly_chart(dna_variant_figure(r["start"], r["end"], variants), use_container_width=True)
        st.caption("DNA 이중나선은 위치를 직관적으로 보여주기 위한 모식도입니다. 실제 분자 구조·나선 회전수는 반영하지 않습니다.")
    else:
        st.success("정렬된 구간에서 염기 차이를 찾지 못했습니다.")

    st.subheader("선택한 ORF의 단백질 영향")
    orfs = st.session_state.get("orfs", [])
    if not orfs:
        st.info("먼저 ③ ORF 분석하기를 눌러 ORF 후보를 생성해주세요.")
        chosen_orf = None
    else:
        recommended = recommend_orfs(orfs, r)
        if recommended:
            st.write("**비교 서열과 겹치는 ORF 추천**")
            st.caption("전체 포함 여부 → ORF 정렬 비율 → 겹치는 염기 수 순으로 추천합니다. "
                       "ORF가 실제 코딩 영역이라는 의미는 아닙니다.")
            best = recommended[:30]
            chosen_index = st.selectbox(
                "단백질 영향 분석에 사용할 ORF",
                [item["index"] for item in best],
                format_func=lambda i: next(
                    f"ORF {i+1} | {item['orf']['start']:,}–{item['orf']['end']:,} bp | "
                    f"길이 {len(item['orf']['dna']):,} bp | 겹침 {item['fraction']:.1%} | "
                    f"{'전체 정렬' if item['full'] else '부분 정렬'}"
                    for item in best if item["index"] == i),
                key="recommended_orf_choice"
            )
            chosen_orf = orfs[chosen_index]
            if not next(item["full"] for item in best if item["index"] == chosen_index):
                st.warning("선택한 ORF가 비교 서열에 완전히 포함되지 않아 단백질 영향 계산이 제한될 수 있습니다.")
        else:
            st.warning("현재 정렬 구간과 겹치는 ORF를 찾지 못했습니다.")
            chosen_orf = None
    effect, reason = analyze_orf_effect(sequence, r, chosen_orf)
    if reason:
        st.info(reason)
    else:
        st.write(effect["summary"])
        if effect["rows"]:
            st.dataframe(pd.DataFrame(effect["rows"]), use_container_width=True, hide_index=True)
        if effect["ref_protein"] != effect["alt_protein"]:
            with st.expander("변이 전·후 예측 아미노산 서열"):
                st.write("기준 아미노산 서열")
                st.code(effect["ref_protein"])
                st.write("비교 아미노산 서열")
                st.code(effect["alt_protein"])
        if effect.get("indels"):
            st.warning("삽입·결실이 포함되어 코돈별 SNP 해석 대신 전체 ORF 번역 결과를 표시합니다. "
                       "프레임시프트 여부는 삽입·결실의 위치와 길이를 추가로 확인해야 합니다.")
    st.caption("이 분석은 선택한 ORF를 코딩 영역으로 가정한 교육용 예측입니다. "
               "실제 단백질 기능 변화나 질병 유발 여부를 판정하지 않습니다. "
               "ClinVar 연동은 참조 서열과 변이 좌표의 정확한 매칭이 필요해 이번 버전에서 제외했습니다.")
