package com.Sts.Collaboration;


public interface SessionListener {

    //// methods
    public void createdSession(SessionManager sm, LeaderChannel lc);

    public void connectedToSession(SessionManager sm, PeerChannel pc);

    public void connectedToNewPeer(SessionManager sm, Participant p);

    public void disconnectedFromSession(SessionManager sm);
    
    public void lostConnectionToPeer(SessionManager sm, Participant p);

    /**
     * This method is called when the SessionManager is promoted to
     * leader of the session.
     */
    public void promoted(SessionManager sm, LeaderChannel lc);

    /**
     * This method is called when the SessionManager is demoted from
     * leader to peer of the session.
     */
    public void demoted(SessionManager sm, PeerChannel pc);
    
    public void peerPromoted(SessionManager sm, Participant p);
}
