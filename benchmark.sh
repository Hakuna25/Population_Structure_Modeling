#!/bin/bash

init_metrics_file() {
    local metrics_file="$1"
    local metrics_dir
    local reset_metrics="${BENCHMARK_RESET_METRICS:-1}"
    metrics_dir="$(dirname "$metrics_file")"
    mkdir -p "$metrics_dir"

    if [[ "$reset_metrics" == "1" || ! -f "$metrics_file" ]]; then
        printf "timestamp\tmethod\tstep\tK\telapsed_sec\tuser_cpu_sec\tsys_cpu_sec\tmax_rss_kb\texit_code\tlog_file\n" > "$metrics_file"
    fi
}

append_metrics_row() {
    local metrics_file="$1"
    local row="$2"
    local lock_dir="${metrics_file}.lockdir"

    while ! mkdir "$lock_dir" 2>/dev/null; do
        sleep 0.1
    done

    printf "%s" "$row" >> "$metrics_file"
    rmdir "$lock_dir"
}

benchmark_run() {
    local metrics_file="$1"
    local method="$2"
    local step="$3"
    local k_value="$4"
    local log_file="$5"
    shift 5

    local start_ts end_ts elapsed_sec exit_code
    local user_cpu_sec="NA"
    local sys_cpu_sec="NA"
    local max_rss_kb="NA"
    local tmp_time
    local timestamp

    start_ts="$(date +%s)"
    tmp_time="$(mktemp)"

    # Clear log file before writing to avoid stale data from previous runs
    : > "$log_file"

    if [[ -x /usr/bin/time ]]; then
        /usr/bin/time -f 'user_cpu_sec=%U\nsys_cpu_sec=%S\nmax_rss_kb=%M' -o "$tmp_time" "$@" \
            > >(tee -a "$log_file") \
            2> >(tee -a "$log_file" >&2)
        exit_code=$?
        user_cpu_sec="$(grep '^user_cpu_sec=' "$tmp_time" | cut -d'=' -f2)"
        sys_cpu_sec="$(grep '^sys_cpu_sec=' "$tmp_time" | cut -d'=' -f2)"
        max_rss_kb="$(grep '^max_rss_kb=' "$tmp_time" | cut -d'=' -f2)"
    else
        "$@" > >(tee -a "$log_file") 2> >(tee -a "$log_file" >&2)
        exit_code=$?
    fi

    end_ts="$(date +%s)"
    elapsed_sec="$((end_ts - start_ts))"
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

    append_metrics_row "$metrics_file" "$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$timestamp" "$method" "$step" "$k_value" "$elapsed_sec" "$user_cpu_sec" "$sys_cpu_sec" "$max_rss_kb" "$exit_code" "$log_file")"

    rm -f "$tmp_time"
    return "$exit_code"
}