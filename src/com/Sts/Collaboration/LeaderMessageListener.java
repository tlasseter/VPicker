package com.Sts.Collaboration;

/**
 * This interface is used by the @link{LeaderChannel} to notify 
 * objects of objects that are received from peers.
 * 
 * @author druths
 */
public interface LeaderMessageListener {

  /**
   * This method is called when an object is received.
   * 
   * @param p is the participant that sent the object.
   * @param obj is the object that was received.
   */
  public void objectReceived(Participant p, Object obj);
}