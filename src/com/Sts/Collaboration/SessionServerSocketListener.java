package com.Sts.Collaboration;

import java.net.*;

/**
 * This interface defines the methods used by the @link{SessionServerSocket} to notify
 * other objects of connection events and status of the SessionServerSocket.
 * 
 * @author druths
 */
public interface SessionServerSocketListener {

    // methods

    // You may need more methods here to handle specific events. But
    // these should give you an idea.

    /**
     * This method is called by the SessionServerSocket when a new connection to
     * it is formed.
     */
    public void connectionFormed(SessionServerSocket sss, Socket s);

    /**
     * This method is called when the SessionServerSocket fails for some reason.
     * This is probably cause to panic the system, since it takes a bad thing to
     * make a server socket die (like loosing network connection).
     */
    public void sessionServerFailed(SessionServerSocket sss);
}