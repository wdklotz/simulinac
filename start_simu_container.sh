function run {
  echo "docker start -a simu"
  docker start -a simu 
}

while true; do
   run
   echo "*: run again"
   echo "q: quit"
   echo "r: reset & quit"
   read RESET
   if [ "$RESET" == "q" ]; then exit; fi
   if [ "$RESET" == "r" ]; then reset; exit; fi
done

