#!/bin/bash
# Author: jwli
# Institution: HZAU
# Date: 2025.9.18

# 显示帮助信息
show_help() {
    cat <<EOF
Usage: $0 [OPTIONS] -p GENEPAIR_LIST -1 CDS1 -2 CDS2 -a PEP1 -b PEP2 -o OUTPUT_DIR

Calculate Ka/Ks ratios for gene pairs using MUSCLE, PAL2NAL, TrimAl/Gblocks and KaKs_Calculator.

Required arguments:
  -p, --genepair FILE    Gene pair list file (tab-delimited)
  -1, --cds1 FILE        CDS sequences for species 1
  -2, --cds2 FILE        CDS sequences for species 2
  -a, --pep1 FILE        Protein sequences for species 1
  -b, --pep2 FILE        Protein sequences for species 2
  -o, --output DIR       Output directory

Optional arguments:
  -h, --help             Show this help message and exit

Example:
  $0 -p genepair.list -1 sample1.cds -2 sample2.cds -a sample1.pep -b sample2.pep -o results
EOF
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_help
            exit 0
            ;;
        -p|--genepair)
            GENEPAIR_LIST="$2"
            shift 2
            ;;
        -1|--cds1)
            SAMPLE1_CDS="$2"
            shift 2
            ;;
        -2|--cds2)
            SAMPLE2_CDS="$2"
            shift 2
            ;;
        -a|--pep1)
            SAMPLE1_PEP="$2"
            shift 2
            ;;
        -b|--pep2)
            SAMPLE2_PEP="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option '$1'"
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$GENEPAIR_LIST" || -z "$SAMPLE1_CDS" || -z "$SAMPLE2_CDS" || \
      -z "$SAMPLE1_PEP" || -z "$SAMPLE2_PEP" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

# 检查软件
command -v muscle >/dev/null || { echo "muscle not found"; exit 1; }
command -v pal2nal.pl >/dev/null || { echo "pal2nal.pl not found"; exit 1; }
command -v trimal >/dev/null || { echo "TrimAl not found"; exit 1; }
command -v Gblocks >/dev/null || { echo "Gblocks not found"; exit 1; }
command -v KaKs_Calculator >/dev/null || { echo "KaKs_Calculator not found"; exit 1; }

# 检查输入文件
for f in "$GENEPAIR_LIST" "$SAMPLE1_CDS" "$SAMPLE2_CDS" "$SAMPLE1_PEP" "$SAMPLE2_PEP"; do
    [ -f "$f" ] || { echo "Error: File $f not found"; exit 1; }
done

mkdir -p "$OUTPUT_DIR"

# 初始化汇总表格文件
SUMMARY_FILE="${OUTPUT_DIR}/kaks_summary.csv"
echo "Gene1,Gene2,Ka,Ks,Ka/Ks,Status" > "$SUMMARY_FILE"

# 初始化合并的kaks结果文件
KAKS_RESULTS="${OUTPUT_DIR}/all_kaks_results.txt"
> "$KAKS_RESULTS"

cleanup() {
    local prefix=$1
    rm -f "${prefix}".{pep,cds,aln,cds.aln,cds.aln.trim,cds.aln.gb,axt,kaks}
}

run_kaks() {
    local axt_file=$1
    local out_file=$2
    local gene1=$3
    local gene2=$4

    KaKs_Calculator -i "$axt_file" -o "$out_file" -m YN

    # 检查KaKs结果是否有效
    if [ -s "$out_file" ] && awk 'NR==2 && $3 != "NA" && $4 != "NA"' "$out_file" | grep -q .; then
        return 0  # 成功
    else
        return 1  # 失败
    fi
}

process_pair() {
    local gene1=$1 gene2=$2
    local cds1=$3 cds2=$4 pep1=$5 pep2=$6
    local outdir=$7
    local tmp_dir="$outdir/tmp"

    mkdir -p "$tmp_dir"
    local prefix="$tmp_dir/${gene1}_${gene2}"

    # 提取序列
    local p1=$(awk -v id="$gene1" -v RS=">" '$1==id {print ">"$0}' "$pep1")
    local p2=$(awk -v id="$gene2" -v RS=">" '$1==id {print ">"$0}' "$pep2")
    local c1=$(awk -v id="$gene1" -v RS=">" '$1==id {print ">"$0}' "$cds1")
    local c2=$(awk -v id="$gene2" -v RS=">" '$1==id {print ">"$0}' "$cds2")

    [ -z "$p1" ] || [ -z "$p2" ] || [ -z "$c1" ] || [ -z "$c2" ] && {
        echo -e "$gene1\t$gene2\tNA\tNA\tNA\tMissing sequence" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
        echo "$gene1,$gene2,NA,NA,NA,Missing sequence" >> "$SUMMARY_FILE"
        return
    }

    echo -e "$p1\n$p2" > "${prefix}.pep"
    echo -e "$c1\n$c2" > "${prefix}.cds"

    # MUSCLE比对
    muscle -align "${prefix}.pep" -output "${prefix}.aln" -threads 2 || {
        echo -e "$gene1\t$gene2\tNA\tNA\tNA\tMUSCLE failed" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
        echo "$gene1,$gene2,NA,NA,NA,MUSCLE failed" >> "$SUMMARY_FILE"
        return
    }

    # PAL2NAL转换
    pal2nal.pl "${prefix}.aln" "${prefix}.cds" -output fasta > "${prefix}.cds.aln" || {
        echo -e "$gene1\t$gene2\tNA\tNA\tNA\tPAL2NAL failed" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
        echo "$gene1,$gene2,NA,NA,NA,PAL2NAL failed" >> "$SUMMARY_FILE"
        return
    }

    # 尝试1: TrimAl处理
    trimal -in "${prefix}.cds.aln" -out "${prefix}.cds.aln.trim" -automated1 || {
        echo -e "$gene1\t$gene2\tNA\tNA\tNA\tTrimAl failed" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
        echo "$gene1,$gene2,NA,NA,NA,TrimAl failed" >> "$SUMMARY_FILE"
        return
    }

    # 转换为AXT格式
    echo "${gene1}_${gene2}" > "${prefix}.axt"
    grep -v "^>" "${prefix}.cds.aln.trim" | grep -v "^$" | grep -v "^-*$" >> "${prefix}.axt"

    # 第一次KaKs计算
    local out="${prefix}.kaks"
    local status="Success"
    local method="TrimAl"

    if ! run_kaks "${prefix}.axt" "$out" "$gene1" "$gene2"; then
        # 第一次计算失败，尝试使用Gblocks

        # 运行Gblocks但不检查其返回值
        Gblocks "${prefix}.cds.aln" -t=c -b4=5 -b5=h -e=.gb

        # 检查Gblocks是否生成了输出文件
        if [ ! -f "${prefix}.cds.aln.gb" ]; then
            echo -e "$gene1\t$gene2\tNA\tNA\tNA\tGblocks output missing" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
            echo "$gene1,$gene2,NA,NA,NA,Gblocks output missing" >> "$SUMMARY_FILE"
            return
        fi

        # 使用Gblocks结果生成新的AXT文件
        echo "${gene1}_${gene2}" > "${prefix}.axt"
        grep -v "^>" "${prefix}.cds.aln.gb" | grep -v "^$" | grep -v "^-*$" |sed 's/ //g' >> "${prefix}.axt"

        # 第二次KaKs计算
        method="Gblocks"
        if ! run_kaks "${prefix}.axt" "$out" "$gene1" "$gene2"; then
            status="KaKs calculation failed after both methods"
            echo -e "$gene1\t$gene2\tNA\tNA\tNA\t$status" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
            echo "$gene1,$gene2,NA,NA,NA,$status" >> "$SUMMARY_FILE"
            return
        else
            status="Success after Gblocks"
        fi
    fi

    # 处理kaks结果文件
    if [ ! -s "$KAKS_RESULTS" ]; then
        # 第一次写入，包含标题
        cat "$out" > "$KAKS_RESULTS"
    else
        # 后续写入，跳过标题
        awk 'NR>1' "$out" >> "$KAKS_RESULTS"
    fi

    local result=$(awk 'NR==2 {print $3"\t"$4"\t"$5}' "$out")
    echo -e "$gene1\t$gene2\t${result}\t$method $status" | tee -a "${OUTPUT_DIR}/kaks_results.txt"
    echo "$gene1,$gene2,$(echo $result | tr '\t' ','),$method $status" >> "$SUMMARY_FILE"

    cleanup "$prefix"
}

# 清空或创建结果文件
> "${OUTPUT_DIR}/kaks_results.txt"

# 主循环
while IFS=$'\t' read -r g1 g2; do
    [[ "$g1" =~ ^# ]] || [[ -z "$g1" || -z "$g2" ]] && continue
    process_pair "$g1" "$g2" "$SAMPLE1_CDS" "$SAMPLE2_CDS" "$SAMPLE1_PEP" "$SAMPLE2_PEP" "$OUTPUT_DIR"
done < "$GENEPAIR_LIST"

echo "Done. Results in:"
echo " - Detailed results: ${OUTPUT_DIR}/kaks_results.txt"
echo " - Summary table: ${SUMMARY_FILE}"
echo " - All KAKS results: ${KAKS_RESULTS}"
