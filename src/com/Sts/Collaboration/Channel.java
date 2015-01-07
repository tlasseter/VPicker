package com.Sts.Collaboration;

/**
 * This interface defines the shared functions between the LeaderChannel and PeerChannel.
 *
 * @author druths
 */
public interface Channel {

    //// Methods
    /**
     * This method returns true only if the channel is active - is
     * managing a valid connection.
     */
    public boolean isActive();
}
