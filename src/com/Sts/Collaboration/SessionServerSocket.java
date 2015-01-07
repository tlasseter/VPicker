package com.Sts.Collaboration;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * This class is a package scope object that allows peers to connect to the
 * local system. It is used by the SessionManager to allow peers to connect.
 * 
 * The SessionManager is notified of peers connecting and errors through the
 * SessionServerSocketListener interface.
 */
class SessionServerSocket {

    // constants
    public static final int DEFAULT_BACKLOG = 50;

    // inner classes
    /**
     * This is the worker thread.  It waits for new connections and notifies the listeners.
     * 
     * @author druths
     */
    private class SSHandler extends Thread {

        public void run() {

            // service the server socket
            try {
                while (true) {

                    Socket s = _ss.accept();

                    notifyListenersOfConnectionFormed(s);
                }
            } catch (IOException ioe) {

                notifyListenersOfFailure();
            }
        }
    }

    // fields
    private HashSet _listeners = new HashSet();

    private ServerSocket _ss;

    private SSHandler _handler;

    // constructors
    public SessionServerSocket(InetAddress addr, int port) throws IOException {
        _ss = new ServerSocket(port, DEFAULT_BACKLOG, addr);

        // start a handling thread
        _handler = new SSHandler();
        _handler.start();
    }

    //// methods
    public void close() {

        _handler.interrupt();

        try {

            _ss.close();
        } catch (IOException ioe) {

            // Don't do anything - we were expecting this.
        }
    }

    // listener methods
    public void addListener(SessionServerSocketListener sssl) {

        _listeners.add(sssl);
    }

    public void removeListener(SessionServerSocketListener sssl) {

        _listeners.remove(sssl);
    }

    protected void notifyListenersOfConnectionFormed(Socket s) {

        Iterator i = _listeners.iterator();

        while (i.hasNext()) {

            SessionServerSocketListener sssl = (SessionServerSocketListener) i
                    .next();

            try {

                sssl.connectionFormed(this, s);

            } catch (RuntimeException rte) {

                System.err.println(sssl
                        + ".connectionFormed threw an exception");
                rte.printStackTrace();

            }
        }
    }

    protected void notifyListenersOfFailure() {

        Iterator i = _listeners.iterator();

        while (i.hasNext()) {

            SessionServerSocketListener sssl = (SessionServerSocketListener) i
                    .next();

            try {

                sssl.sessionServerFailed(this);

            } catch (RuntimeException rte) {

                System.err.println(sssl
                        + ".sessionServerFailed threw an exception");
                rte.printStackTrace();

            }
        }
    }
}