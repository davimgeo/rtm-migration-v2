setsid make &
MAKE_PID=$!

sleep 1

kill -- -$MAKE_PID
wait $MAKE_PID 2>/dev/null

bash send_to_thesis.sh

cd ~/coding/undergraduate-thesis/
make
