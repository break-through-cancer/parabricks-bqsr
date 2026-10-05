def pbrunFunction() {
    '''
pbrun() {
    local rc=0 exe
    exe=$(type -P pbrun)
    "$exe" "$@" 2>&1 | tee -a pbrun.log || rc=$?
    if [ "$rc" -ne 0 ] && grep -q 'SIGKILL' pbrun.log; then
        echo "pbrun was killed for host memory (SIGKILL); exiting 137 so the task is retried with more memory" >&2
        return 137
    fi
    return "$rc"
}
'''
}
