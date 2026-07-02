#!/bin/bash
# Nights/weekends driver for the d=5 campaign.  Runs the resumable queue with:
#   - caffeinate -i : prevent IDLE sleep while working (still sleeps on lid close)
#   - nice          : lowest priority, never fights an interactive session
#   - a wall-clock deadline: weekdays -> stop at STOP_HOUR (default 07:00); weekends -> run long
# Safe to launch repeatedly (launchd/cron or by hand); the filesystem queue resumes.
#
#   ./run_d5_campaign_night.sh            # auto: until 7am on weekdays, 14h on weekends
#   ./run_d5_campaign_night.sh 120        # explicit: run 120 minutes
set -euo pipefail
cd "$(dirname "$0")"

NPROC="${NPROC:-7}"
STOP_HOUR="${STOP_HOUR:-7}"

if [ "${1:-}" != "" ]; then
  MINUTES="$1"
else
  DOW=$(date +%u)          # 1=Mon .. 7=Sun
  if [ "$DOW" -ge 6 ]; then
    MINUTES=840            # weekend: 14h chunk (relaunch continues)
  else
    now_h=$(date +%H); now_m=$(date +%M)
    # minutes from now until STOP_HOUR (next occurrence)
    MINUTES=$(( ( (STOP_HOUR*60) - (10#$now_h*60 + 10#$now_m) + 1440 ) % 1440 ))
    [ "$MINUTES" -lt 15 ] && MINUTES=$((MINUTES+1440))   # if already past, go to tomorrow
  fi
fi

echo "$(date '+%F %T')  campaign session: ${MINUTES} min, ${NPROC} workers"
exec caffeinate -i nice -n 15 python3 run_d5_campaign.py run "$MINUTES" "$NPROC"
