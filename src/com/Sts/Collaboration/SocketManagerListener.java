package com.Sts.Collaboration;

/**
 * This interface defines the methods used by the @link{SocketManager} to notify other objects of
 * it state.
 * 
 * @author druths
 */
public interface SocketManagerListener {
    
    public void socketManagerDeactivated(SocketManager sm);
}
