import streamlit as st
import pandas as pd
import json
import os
import time
import subprocess
import sys
import html as html_lib
from datetime import datetime
from core_pipeline import run_pipeline
from pyhpo import Ontology

# ============================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================
def generate_population_report(results_dict):
    lines = ["gnomAD v4 (AF вычисляется как ac/an)\n"]
    df = results_dict.get('primary')
    if df is None or df.empty:
        lines.append("Убедительных кандидатов не найдено.")
        return "\n".join(lines)
    for _, row in df.iterrows():
        gene, variant, af, score = row['gene'], row['variant'], row.get('af'), row['score']
        if score == 0: continue
        lines.append(f"--- {gene} ({variant}) ---")
        if af is None:
            lines.append(f"  не найден в gnomAD (крайне редкий).")
        elif af > 0.01:
            lines.append(f"  Частый вариант (AF = {af:.6f}), вероятно доброкачественный.")
        else:
            lines.append(f"  Редкий вариант (AF = {af:.6f}), может быть патогенным.")
        lines.append("")
    return "\n".join(lines)

def generate_html_report(results_dict, filename):
    now = datetime.now().strftime("%d.%m.%Y %H:%M")
    def esc(x): return html_lib.escape(str(x)) if x is not None else ""
    
    html = f"""<html><head><meta charset="utf-8"><title>Genomic Detective Report</title>
    <style>body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; line-height: 1.6; }}
    h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 10px; }}
    h2 {{ color: #2c3e50; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 10px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 15px; }}
    th, td {{ border: 1px solid #ddd; padding: 10px; text-align: left; }}
    th {{ background-color: #f8f9fa; }}
    .pathogenic {{ color: #c0392b; font-weight: bold; background-color: #fadbd8; padding: 2px 6px; border-radius: 4px; }}
    .recommendation {{ background-color: #e8f4f8; padding: 20px; border-left: 5px solid #3498db; margin-top: 20px; }}
    .disclaimer {{ background-color: #fff3cd; color: #856404; padding: 15px; border-left: 5px solid #ffeeba; margin-bottom: 20px; }}
    .footer {{ margin-top: 50px; font-size: 12px; color: #7f8c8d; text-align: center; border-top: 1px solid #eee; padding-top: 20px; }}</style></head><body>
    <h1>🧬 Отчет Genomic Detective</h1>
    <p><strong>Дата:</strong> {now}<br><strong>Файл:</strong> {esc(filename)}</p>
    <div class="disclaimer"><strong>⚠️ ВНИМАНИЕ:</strong> Требуется валидация в аккредитованной лаборатории.</div>"""
    
    df = results_dict.get('primary')
    if df is not None and not df.empty:
        html += "<h2>🩺 Приоритетные гены-кандидаты</h2><table><tr><th>Ранг</th><th>Ген</th><th>Мутация</th><th>Болезнь</th><th>ClinVar</th><th>Балл</th></tr>"
        for i, r in df.head(5).iterrows():
            cv_class = "pathogenic" if r['clinvar'] in ("pathogenic", "likely_pathogenic") else ""
            html += f"<tr><td>{i+1}</td><td><b>{esc(r['gene'])}</b></td><td>{esc(r['variant'])}</td><td>{esc(r['disease'])}</td><td class='{cv_class}'>{esc(r['clinvar'])}</td><td>{r['score']}</td></tr>"
        html += "</table>"
        
        top1 = df.iloc[0]
        html += f"<div class='recommendation'><h3>Рекомендация по кандидату {esc(top1['gene'])}</h3>"
        html += f"<p><b>Тип наследования:</b> {esc(top1['inheritance'])}<br><b>Балл:</b> {top1['score']}</p>"
        html += f"<p>{esc(top1.get('recommendation', 'Рекомендована валидация.'))}</p></div>"
        
    html += "<div class='footer'>Сгенерировано Genomic Detective. Не является медицинским заключением.</div></body></html>"
    return html

# ============================================================
# НАСТРОЙКА СТРАНИЦЫ
# ============================================================
st.set_page_config(page_title="Genomic Detective", page_icon="🧬", layout="wide")
st.warning("⚠️ **ДЕМОНСТРАЦИОННАЯ ВЕРСИЯ.** Результаты носят исключительно исследовательский характер. Не является клиническим диагнозом.")
st.title("🧬 Genomic Detective 🧬")
st.markdown("Загрузите данные пациента для проведения анализа")

@st.cache_resource
def load_hpo_terms():
    Ontology()
    terms = {}
    for term in Ontology:
        if "inheritance" in term.name.lower(): continue
        if term.id not in ["HP:0000001", "HP:0000118", "HP:0000005"]:
            terms[f"{term.name} ({term.id})"] = term.id
    return terms

hpo_dict = load_hpo_terms()

if "saved_reports" not in st.session_state: st.session_state.saved_reports = []
if "current_results" not in st.session_state: st.session_state.current_results = None
if "current_filename" not in st.session_state: st.session_state.current_filename = ""

# ============================================================
# БОКОВАЯ ПАНЕЛЬ
# ============================================================
with st.sidebar:
    st.header("⚙️ Настройки")
    st.header("1. Данные о вариантах")
    data_format = st.radio("Формат данных:", ["Готовая таблица (TSV)", "SNP (23andMe / AncestryDNA)", "Сырые риды (FASTQ - WGS)"])

    variant_file = None
    snp_file = None
    r1_file, r2_file = None, None

    if data_format == "Готовая таблица (TSV)":
        variant_file = st.file_uploader("Файл вариантов (.tsv)", type=['tsv'])
    elif data_format == "SNP (23andMe / AncestryDNA)":
        snp_file = st.file_uploader("Файл SNP (.txt / .csv)", type=['txt', 'csv'])
    else:
        r1_file = st.file_uploader("Файл R1 (.fastq.gz)", type=['gz', 'fastq'])
        r2_file = st.file_uploader("Файл R2 (.fastq.gz)", type=['gz', 'fastq'])

    st.header("2. Фенотип пациента")
    pheno_mode = st.radio("Способ ввода симптомов:", ["Заполнить вручную", "Загрузить файл с симптомами"])

    json_pheno_file = None
    selected_symptoms = []
    excluded_symptoms = []

    if pheno_mode == "Загрузить файл с симптомами":
        json_pheno_file = st.file_uploader("Файл симптомов (.json)", type=['json'])
    else:
        selected_symptoms = st.multiselect("Подтвержденные симптомы:", options=list(hpo_dict.keys()))
        excluded_symptoms = st.multiselect("Исключенные симптомы:", options=list(hpo_dict.keys()))

    st.header("3. Происхождение пациента")
    meta_mode = st.radio("Способ указания популяции:", ["Выбрать вручную", "Загрузить файл метаданных"])
    population = "nfe"
    json_meta_file = None

    if meta_mode == "Выбрать вручную":
        population = st.selectbox("Популяция (для gnomAD):", ["nfe", "afr", "eas", "sas", "amr", "asj", "fin", "mid"])
    else:
        json_meta_file = st.file_uploader("Файл метаданных (.json)", type=['json'], help="Ожидается поле 'ancestry'")
        if json_meta_file:
            try:
                meta_data = json.load(json_meta_file)
                json_meta_file.seek(0)
                detected_pop = meta_data.get("ancestry", meta_data.get("population", "nfe"))
                st.success(f"Определена популяция: {detected_pop}")
                population = detected_pop
            except Exception as e:
                st.error(f"Ошибка чтения метаданных: {e}")

    st.markdown("---")
    run_btn = st.button("🔍 Запустить анализ", type="primary", use_container_width=True)

# ============================================================
# ОСНОВНОЙ ИНТЕРФЕЙС
# ============================================================
tab1, tab2, tab3, tab4 = st.tabs(["🩺 Основной анализ", "🧬 Носительство", "⚠️ Вторичные находки", "📊 Сравнение"])

with tab1:
    if run_btn:
        if data_format == "Готовая таблица (TSV)" and not variant_file:
            st.error("⚠️ Загрузите файл вариантов (.tsv).")
        elif data_format == "SNP (23andMe / AncestryDNA)" and not snp_file:
            st.error("⚠️ Загрузите файл SNP (.txt / .csv).")
        elif data_format == "Сырые риды (FASTQ - WGS)" and (not r1_file or not r2_file):
            st.error("⚠️ Загрузите оба файла R1 и R2.")
        elif pheno_mode == "Загрузить файл с симптомами" and not json_pheno_file:
            st.error("⚠️ Загрузите файл симптомов (.json).")
        elif pheno_mode == "Заполнить вручную" and not selected_symptoms:
            st.error("⚠️ Выберите хотя бы один симптом.")
        else:
            temp_variant_path, temp_pheno_path = "temp_variants.tsv", "temp_pheno.json"
            try:
                # Сохранение фенотипа
                if json_pheno_file:
                    with open(temp_pheno_path, "wb") as f: f.write(json_pheno_file.getbuffer())
                else:
                    features = [{"type": {"id": hpo_dict[s]}, "excluded": False} for s in selected_symptoms] + \
                               [{"type": {"id": hpo_dict[s]}, "excluded": True} for s in excluded_symptoms]
                    with open(temp_pheno_path, "w", encoding='utf-8') as f: json.dump({"phenotypicFeatures": features}, f)

                if data_format == "Готовая таблица (TSV)":
                    with open(temp_variant_path, "wb") as f: f.write(variant_file.getbuffer())
                    with st.spinner("Валидация и расчет HPO-сходства..."):
                        st.session_state.current_results = run_pipeline(temp_variant_path, temp_pheno_path, population=population)
                        st.session_state.current_filename = variant_file.name

                else:
                    # WGS или SNP через subprocess
                    workdir = "work/current_run"
                    os.makedirs(workdir, exist_ok=True)

                    if data_format == "SNP (23andMe / AncestryDNA)":
                        file1 = os.path.join(workdir, "input.txt")
                        with open(file1, "wb") as f: f.write(snp_file.getbuffer())
                        type_arg = "snp"
                        cmd = [sys.executable, "run_pipeline.py", type_arg, file1, workdir]
                    else:
                        file1 = os.path.join(workdir, "input1.fastq.gz")
                        file2 = os.path.join(workdir, "input2.fastq.gz")
                        with open(file1, "wb") as f: f.write(r1_file.getbuffer())
                        with open(file2, "wb") as f: f.write(r2_file.getbuffer())
                        type_arg = "wgs"
                        cmd = [sys.executable, "run_pipeline.py", type_arg, file1, file2, workdir]

                    with open(os.path.join(workdir, "pheno.json"), "w") as f:
                        f.write(open(temp_pheno_path).read())

                    proc = subprocess.Popen(cmd)
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    while proc.poll() is None:
                        status_file = os.path.join(workdir, "status.json")
                        if os.path.exists(status_file):
                            with open(status_file) as f:
                                status = json.load(f)
                            progress_bar.progress(status.get("progress", 0) / 100)
                            status_text.info(f"**{status.get('step', '')}**: {status.get('message', '')}")
                        time.sleep(1)

                    if proc.returncode == 0:
                        result_file = os.path.join(workdir, "result.json")
                        if os.path.exists(result_file):
                            with open(result_file) as f:
                                res_data = json.load(f)
                            df_primary = pd.DataFrame(res_data['ranked'])
                            st.session_state.current_results = {"primary": df_primary, "carrier": pd.DataFrame(), "secondary": pd.DataFrame()}
                            st.session_state.current_filename = "WGS_processed" if type_arg == "wgs" else snp_file.name
                            st.success("✅ Анализ завершен!")
                        else:
                            st.error("❌ Результат не найден.")
                    else:
                        st.error("❌ Ошибка при выполнении пайплайна. Проверьте логи.")
                        st.stop()
            except Exception as e:
                st.error(f"❌ Ошибка: {e}")
            finally:
                if os.path.exists(temp_variant_path): os.remove(temp_variant_path)
                if os.path.exists(temp_pheno_path): os.remove(temp_pheno_path)

    # --- ОТОБРАЖЕНИЕ РЕЗУЛЬТАТОВ С КРАСИВЫМИ КАРТОЧКАМИ ---
    if st.session_state.current_results is not None:
        results_dict = st.session_state.current_results
        if results_dict['primary'].empty:
            st.error("Убедительных кандидатов не найдено.")
        else:
            st.success(f"Анализ для файла **{st.session_state.current_filename}** завершен!")
            st.subheader("Приоритетные гены-кандидаты")

            for i in range(min(5, len(results_dict['primary']))):
                row = results_dict['primary'].iloc[i]
                status_color = "🔴" if "pathogenic" in str(row['clinvar']).lower() else "🟡" if "uncertain" in str(row['clinvar']).lower() else "⚪"
                header = f"{status_color} {i+1}. Ген: {row['gene']} | Мутация: {row['variant']} | Балл (log10 LR): {row['score']}"

                with st.expander(header, expanded=(i==0)):
                    # 🚨 ВИЗУАЛЬНЫЙ АПГРЕЙД
                    c1, c2, c3 = st.columns(3)
                    c1.metric("Итоговый балл", row['score'])
                    c2.metric("Совпадение с HPO", row.get('hpo_sim', 'N/A'))
                    penalty_val = row.get('penalty', 0)
                    c3.metric("Штраф за исключенные", f"-{penalty_val}", delta=f"-{penalty_val}", delta_color="inverse")
                    
                    st.markdown("---")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("### 🧬 Данные о варианте")
                        st.write(f"**Тип мутации:** {row['consequence']}")
                        st.write(f"**Статус в ClinVar:** {row['clinvar']}")
                        st.write(f"**Тип наследования:** {row['inheritance']}")
                    with col2:
                        st.markdown("### 🩺 Клиническая картина")
                        st.write(f"**Связанное заболевание:** {row['disease']}")
                        if 'af' in row and pd.notna(row['af']):
                            st.write(f"**Частота в gnomAD:** {row['af']}")
                        else:
                            st.write(f"**Частота в gnomAD:** не найден (редкий)")
                    
                    st.info(f"**Рекомендация:** {row.get('recommendation', 'Рекомендована валидация методом Сэнгера.')}")

            st.markdown("---")
            st.subheader("📊 Анализ частот в gnomAD (текстовый отчет)")
            report_text = generate_population_report(results_dict)
            with st.expander("Показать текстовый отчет", expanded=False):
                st.text(report_text)

            st.markdown("---")
            st.markdown("### 📥 Экспорт результатов")
            col_dl1, col_dl2, col_dl3 = st.columns(3)
            
            # HTML
            html_report = generate_html_report(results_dict, st.session_state.current_filename)
            col_dl1.download_button("📄 Клинический отчет (HTML)", html_report, file_name=f"report_{st.session_state.current_filename}.html", mime="text/html", use_container_width=True)
            # CSV
            col_dl2.download_button("📊 Сырые данные (CSV)", data=results_dict['primary'].to_csv(index=False).encode('utf-8'), file_name=f"report_{st.session_state.current_filename}.csv", mime="text/csv", use_container_width=True)
            # TXT
            col_dl3.download_button("📥 Частоты (TXT)", data=report_text, file_name="population_analysis.txt", mime="text/plain", use_container_width=True)

            if st.button("💾 Сохранить в историю для сравнения"):
                st.session_state.saved_reports.append({"filename": st.session_state.current_filename, "df": results_dict['primary'].copy()})
                st.toast("✅ Отчет сохранен! Перейдите во вкладку 'Сравнение'.")
    else:
        st.info("👈 Заполните данные в левой панели и нажмите «Запустить анализ».")

with tab2:
    st.subheader("Анализ носительства")
    if st.session_state.current_results is not None and not st.session_state.current_results['carrier'].empty:
        st.warning("Выявлены варианты, которые могут передаться потомству.")
        st.dataframe(st.session_state.current_results['carrier'], use_container_width=True, hide_index=True)
    else:
        st.info("Носительство не обнаружено или анализ не запущен.")

with tab3:
    st.subheader("Вторичные находки (ACMG SF)")
    if st.session_state.current_results is not None and not st.session_state.current_results['secondary'].empty:
        st.error("Обнаружены клинически значимые варианты.")
        st.dataframe(st.session_state.current_results['secondary'], use_container_width=True, hide_index=True)
    else:
        st.info("Вторичных находок не обнаружено или анализ не запущен.")

with tab4:
    st.subheader("Сравнение сохраненных отчетов")
    if "report_to_delete" not in st.session_state: st.session_state.report_to_delete = None

    if not st.session_state.saved_reports:
        st.info("Нет сохраненных отчетов. Проведите анализ и нажмите «Сохранить в историю для сравнения».")
    else:
        if st.session_state.report_to_delete is not None:
            del st.session_state.saved_reports[st.session_state.report_to_delete]
            st.session_state.report_to_delete = None
            st.rerun()

        comparison_data = {}
        for report in st.session_state.saved_reports:
            comparison_data[report['filename']] = report['df'].set_index('gene')['score']

        if comparison_data:
            comparison_df = pd.DataFrame(comparison_data).sort_index().fillna("-")
            st.markdown("**Сводная таблица баллов по генам (log10 LR)**")
            st.dataframe(comparison_df, use_container_width=True)

            st.markdown("---")
            st.markdown("**Детали по отчетам**")

            for idx, report in enumerate(st.session_state.saved_reports):
                col1, col2 = st.columns([0.9, 0.1])
                with col1:
                    with st.expander(f"📄 Отчет {idx+1}: {report['filename']}"):
                        st.dataframe(report['df'][['gene', 'variant', 'score', 'clinvar']].head(5), use_container_width=True, hide_index=True)
                with col2:
                    if st.button("🗑️", key=f"del_{idx}", help="Удалить этот отчет"):
                        st.session_state.report_to_delete = idx
                        st.rerun()