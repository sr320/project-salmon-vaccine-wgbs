#!/usr/bin/env bash
# Cron trigger: start the M-bias pass (code/04-dedup-methylation.sh mbias)
# once Bismark alignment (code/03-bismark.sh) has finished.
#
# Installed as a cron entry tagged '# salmon-wgbs-mbias-trigger' that runs
# every 15 minutes. On each run it:
#   - does nothing while alignment is still going;
#   - if alignment finished, removes its own cron entry and runs the M-bias pass;
#   - if alignment stopped without finishing, removes its cron entry and logs why.
set -uo pipefail

REPO=/home/shared/8TB_HDD_03/sr320/github/project-salmon-vaccine-wgbs
ALN_LOG=$REPO/output/03-bismark/03-bismark.log
LOG=$REPO/output/04-methylation/04-mbias-trigger.log
TAG='# salmon-wgbs-mbias-trigger'
export PATH=/usr/local/bin:/usr/bin:/bin

remove_cron() { ( crontab -l 2>/dev/null | grep -v "$TAG" ) | crontab -; }

# Only one copy at a time.
exec 9>"/tmp/salmon-wgbs-mbias-trigger.lock"
flock -n 9 || exit 0

if grep -q "== finished" "$ALN_LOG" 2>/dev/null; then
  mkdir -p "$(dirname "$LOG")"
  remove_cron
  echo "$(date): alignment finished; starting M-bias pass" >> "$LOG"
  "$REPO/code/04-dedup-methylation.sh" mbias >> "$LOG" 2>&1
  echo "$(date): M-bias pass exited with status $?" >> "$LOG"
elif ! pgrep -f "code/03-bismark.sh" >/dev/null; then
  mkdir -p "$(dirname "$LOG")"
  remove_cron
  echo "$(date): alignment job is not running and did not finish; M-bias not started" >> "$LOG"
fi
