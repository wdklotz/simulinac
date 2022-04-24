import zmq
import time
import sys
import random

port = "5560"
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.connect("tcp://localhost:%s" % port)
server_id = random.randrange(1,10005)
while True:
    #  Wait for next request from client
    message = socket.recv()
    print("Received request: {}".format(message))
    time.sleep (0.1)  
    socket.send_string("<WELCOME> from server %s" % server_id)
