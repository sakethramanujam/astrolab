#!/usr/bin/env bash
# One thermal-guard tick (run every minute by a user timer, launchd, or a loop). Reads the CPU temperature without
# sudo. At or above HOT it pauses (SIGSTOP) the astrolab v1.1 workers; below COOL it resumes (SIGCONT) the ones it
# paused. It never kills anything and never touches other processes. State and the last reading go to STATUS (JSON).
#   thermal_guard.sh STATUS_FILE        env: HOT=90 COOL=80 PATTERN='exosearch.v11.run work'
# Linux: hwmon coretemp/k10temp/zenpower, else thermal_zone x86_pkg_temp. macOS: `pmset -g therm` (thermal or
# performance warning recorded, or CPU_Speed_Limit < 80, counts as serious). WSL: powershell.exe
# MSAcpi_ThermalZoneTemperature if Windows allows it unelevated; otherwise no reading and no action.
set -uo pipefail
STATUS="$1"; HOT="${HOT:-90}"; COOL="${COOL:-80}"; PATTERN="${PATTERN:-exosearch.v11.run work}"
STATE="$STATUS.paused"; mkdir -p "$(dirname "$STATUS")"
src=""; temp=""; level=""
if [ "$(uname)" = Darwin ]; then
  src="pmset"; t=$(pmset -g therm 2>/dev/null)
  lim=$(echo "$t" | awk -F= '/CPU_Speed_Limit/ {gsub(/ /,"",$2); print $2}')
  if echo "$t" | grep -q -E "^(Thermal|Performance) warning level" || { [ -n "$lim" ] && [ "$lim" -lt 80 ]; }; then level=serious
  else level=nominal; fi
elif grep -qi microsoft /proc/version 2>/dev/null; then
  src="MSAcpi_ThermalZoneTemperature"
  temp=$(timeout 20 powershell.exe -NoProfile -Command '$t = Get-CimInstance -Namespace root/wmi -ClassName MSAcpi_ThermalZoneTemperature -ErrorAction Stop; [math]::Round((($t | Measure-Object CurrentTemperature -Maximum).Maximum / 10 - 273.15), 1)' 2>/dev/null | tr -d '\r' | grep -E '^[0-9.]+$')
else
  for h in /sys/class/hwmon/hwmon*; do
    case "$(cat "$h/name" 2>/dev/null)" in coretemp|k10temp|zenpower)
      src="hwmon:$(cat "$h/name")"; temp=$(cat "$h"/temp*_input 2>/dev/null | sort -n | tail -1);; esac
  done
  if [ -z "$temp" ]; then
    for z in /sys/class/thermal/thermal_zone*; do
      [ "$(cat "$z/type" 2>/dev/null)" = x86_pkg_temp ] && { src="x86_pkg_temp"; temp=$(cat "$z/temp"); }
    done
  fi
  [ -n "$temp" ] && temp=$((temp / 1000))
fi
if [ -n "$temp" ]; then
  hot=$(awk -v t="$temp" -v h="$HOT" 'BEGIN{print (t>=h)}'); cool=$(awk -v t="$temp" -v c="$COOL" 'BEGIN{print (t<c)}')
elif [ -n "$level" ]; then
  hot=$([ "$level" = serious ] && echo 1 || echo 0); cool=$([ "$level" = nominal ] && echo 1 || echo 0)
else
  hot=0; cool=0                                             # no reading: never act
fi
pids=""; for p in $(pgrep -f "$PATTERN"); do pids="$pids $p $(pgrep -P "$p" | tr '\n' ' ')"; done
action=none
if [ "$hot" = 1 ] && [ -n "${pids// /}" ]; then
  kill -STOP $pids 2>/dev/null; echo "$pids" > "$STATE"; action=paused
elif [ "$cool" = 1 ] && [ -f "$STATE" ]; then
  kill -CONT $(cat "$STATE") $pids 2>/dev/null; rm -f "$STATE"; action=resumed
elif [ -f "$STATE" ]; then
  action=held
fi
printf '{"at_utc": "%s", "source": "%s", "temp_c": %s, "level": "%s", "hot_c": %s, "cool_c": %s, "workers": %s, "action": "%s", "paused": %s}\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${src:-none}" "${temp:-null}" "${level:-}" "$HOT" "$COOL" \
  "$(echo $pids | wc -w | tr -d ' ')" "$action" "$([ -f "$STATE" ] && echo true || echo false)" > "$STATUS"
[ "$action" = paused ] || [ "$action" = resumed ] && cat "$STATUS" >> "$STATUS.log"
exit 0
