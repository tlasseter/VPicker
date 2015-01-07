/*
 * Created on May 16, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package com.Sts.Collaboration;

import java.io.*;

/**
 * This interface is used by the @link{PeerChannel} in order to report new messages received
 * from the leader.
 * 
 * @author druths
 */
public interface PeerMessageListener {
    
    /**
     * This method is called each time a message is received by the @link{PeerChannel}.
     * 
     * @param obj is the object that was received.
     */
    public void messageReceived(Serializable obj);
}
