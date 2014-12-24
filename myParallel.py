#!/usr/bin/python

import Queue
import threading
import time
import sys

exitFlag = 0
queueLock = threading.Lock()


class doThread (threading.Thread):
    def __init__(self, threadID, q, FUN=None, *args):
        threading.Thread.__init__(self)
        self.threadID = threadID
        #self.queue = q
        self.FUN = FUN
        self.args = args
        print >>sys.stderr, args
        if FUN is None:
            print >> sys.stderr, "cannot run a thread with a null function"
            sys.exit(1)

    def run(self):
        print "Starting thread " + self.name
        while not (exitFlag or workQueue.empty()):
        #while not exitFlag:
            queueLock.acquire() #lock Q, ensure data only taken once
            if workQueue.empty():
                queueLock.release()
                print "Exiting " + self.name
            else:
#                data = self.queue.get()
                data= workQueue.get()
                queueLock.release()
                print "%s processing %s\n" % (self.threadID, data)
                results = self.FUN(data, self.args)
            #sleep(5) #if queue empty, wait a second for others to catch up
                
class writeThread (threading.Thread):
    def __init__(self, outfile):
        threading.Thread.__init__(self)
        # file writing daemon - not needed right now, but write one.
        pass
    def run(self):
        pass

def parallel(FUN, noThreads, workArray, *args):
    #fill queue with blocks
    global queueLock
    global workQueue
    # global resultsQueue
    global exitFlag
    
    workQueue = Queue.Queue(len(workArray))
    resultqueue = Queue.Queue()
    queueLock.acquire()
    for block in workArray:
        workQueue.put(block)
    queueLock.release()

    # get threads array - will need it to kill later...
    threads = []
    
    # Create new threads (one less than total threads)
    for tID in range(1,noThreads):
        thread = doThread(tID, workQueue, FUN, *args)
        thread.start()
        threads.append(thread)

    # Wait for queue to empty
    while not workQueue.empty():
        #  make current thread into write-thread
        pass
    # kill other threads if head detects empty queue
    exitFlag = 1
    
    # Wait for all threads to complete
    for t in threads:
        t.join()
    print "Exiting Parallel"
