package com.Sts.Collaboration;

import junit.framework.*;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * @author druths
 */
public class SessionManagerConnectionTest extends TestCase
{

	private class TestMsgHandler implements PeerMessageListener
	{

		public LinkedList queue = new LinkedList();

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.MessageListener#messageReceived(java.io.Serializable)
		 */
		public void messageReceived(Serializable obj)
		{
			queue.addLast(obj);
		}
	}

	private class TestSessionListener implements SessionListener
	{

		// fields
		public boolean connected = false;

		public PeerChannel pc;

		public boolean created = false;

		public LeaderChannel lc;

		public boolean promoted = false;

		public boolean demoted = false;

		public Vector participants = new Vector();

		// methods
		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#connectedToNewPeer(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.Participant)
		 */
		public void connectedToNewPeer(SessionManager sm, Participant p)
		{
			participants.add(p);
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#connectedToSession(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.PeerChannel)
		 */
		public void connectedToSession(SessionManager sm, PeerChannel pc)
		{
			connected = true;
			this.pc = pc;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#createdSession(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.LeaderChannel)
		 */
		public void createdSession(SessionManager sm, LeaderChannel lc)
		{
			created = true;
			this.lc = lc;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#demoted(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.PeerChannel)
		 */
		public void demoted(SessionManager sm, PeerChannel pc)
		{
			demoted = true;
			this.pc = pc;
			System.out.println(sm.getIdentity() + " demoted.");
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#disconnectedFromSession(com.s2s.sm.SessionManager)
		 */
		public void disconnectedFromSession(SessionManager sm)
		{
			// TODO Auto-generated method stub
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#promoted(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.LeaderChannel)
		 */
		public void promoted(SessionManager sm, LeaderChannel lc)
		{
			promoted = true;
			this.lc = lc;
			System.out.println(sm.getIdentity() + " promoted.");
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#lostConnectionToPeer(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.Participant)
		 */
		public void lostConnectionToPeer(SessionManager sm, Participant p)
		{
			// TODO Auto-generated method stub
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionListener#peerPromoted(com.s2s.sm.SessionManager,
		 *      com.s2s.sm.Participant)
		 */
		public void peerPromoted(SessionManager sm, Participant p)
		{
			// TODO Auto-generated method stub
		}
	}

	public void testPeerLeaderPromotion() throws Exception
	{

		// create a leader
		SessionManager leader = new SessionManager();
		SessionManager peer = new SessionManager();
		SessionManager tmp_sm;

		// set listeners
		TestSessionListener peer_listener = new TestSessionListener();
		TestSessionListener leader_listener = new TestSessionListener();
		TestSessionListener tmp_listener;

		leader.addSessionListener(leader_listener);
		peer.addSessionListener(peer_listener);

		// connect them
		leader.createSession(InetAddress.getLocalHost(), 9010, "password");
		peer.connectToSession(InetAddress.getLocalHost(), 9010, 9011, "password");

		// check listeners
//		Assert.assertTrue(leader_listener.created);
//		Assert.assertTrue(peer_listener.connected);

		// do a promotion
		leader_listener.lc.demote(peer.getIdentity());

		Thread.sleep(10000);

//		Assert.assertTrue(leader_listener.demoted);
//		Assert.assertTrue(peer_listener.promoted);

		// swap roles
		tmp_sm = leader;
		leader = peer;
		peer = tmp_sm;

		tmp_listener = leader_listener;
		leader_listener = peer_listener;
		peer_listener = tmp_listener;

		// pass some messages
		TestMsgHandler ml = new TestMsgHandler();
		peer_listener.pc.addListener(ml);

		for (int i = 0; i < 10; i++)
		{

			leader_listener.lc.broadcast(null, new Integer(i));
		}

		Thread.sleep(1000);

		System.out.println("OUT SIZE = " + ml.queue.size());

		int prev = ( (Integer) ml.queue.removeFirst()).intValue();
		int curr;
		while (ml.queue.size() == 0)
		{

			curr = ( (Integer) ml.queue.removeFirst()).intValue();
			System.out.println("INT = " + curr + " PREV = " + prev);
			Assert.assertEquals(curr, prev + 1);

			prev = curr;
		}

		// disconnect
		peer.disconnect();

		Thread.sleep(1000);

		try
		{

			leader_listener.lc.broadcast(null, new Vector());
		}
		catch (Exception e)
		{

			e.printStackTrace();
		}

		leader.disconnect();
	}

	public void testPeerLeaderConnection() throws Exception
	{

		// create a leader
		SessionManager leader = new SessionManager();
		SessionManager peer = new SessionManager();

		// set listeners
		TestSessionListener peer_listener = new TestSessionListener();
		TestSessionListener leader_listener = new TestSessionListener();

		leader.addSessionListener(leader_listener);
		peer.addSessionListener(peer_listener);

		// connect them
		leader.createSession(InetAddress.getLocalHost(), 9010, "password");
		peer.connectToSession(InetAddress.getLocalHost(), 9010, 9011, "password");

		// check listeners
		Assert.assertTrue(leader_listener.created);
		Assert.assertTrue(peer_listener.connected);

		// pass some messages
		TestMsgHandler ml = new TestMsgHandler();
		peer_listener.pc.addListener(ml);

		for (int i = 0; i < 10; i++)
		{

			leader_listener.lc.broadcast(null, new Integer(i));
		}

		Thread.sleep(1000);

		System.out.println("OUT SIZE = " + ml.queue.size());

		int prev = ( (Integer) ml.queue.removeFirst()).intValue();
		int curr;
		while (ml.queue.size() == 0)
		{

			curr = ( (Integer) ml.queue.removeFirst()).intValue();
			System.out.println("INT = " + curr + " PREV = " + prev);
			Assert.assertEquals(curr, prev + 1);

			prev = curr;
		}

		// disconnect
		peer.disconnect();

		Thread.sleep(1000);

		try
		{

			leader_listener.lc.broadcast(null, new Vector());
		}
		catch (Exception e)
		{

			e.printStackTrace();
		}

		leader.disconnect();
	}

	public void testNPeers() throws Exception
	{

		// create one leader and N participants
		int N = 5;
		int lp = 9011;

		Thread.sleep(2000);

		SessionManager leader = new SessionManager();
		SessionManager[] peers = new SessionManager[N];
		TestSessionListener[] tsls = new TestSessionListener[N + 1];

		// start session
		tsls[N] = new TestSessionListener();
		leader.addSessionListener(tsls[N]);
		leader.createSession(InetAddress.getLocalHost(), 9010, "pwd");

		for (int i = 0; i < N; i++)
		{
			System.out.println("********* i = " + i);
			peers[i] = new SessionManager();
			tsls[i] = new TestSessionListener();

			peers[i].addSessionListener(tsls[i]);
			peers[i].connectToSession(InetAddress.getLocalHost(), 9010, lp++, "pwd");

			// let the session stablize
			Thread.sleep(1000);
		}

		// verify that all SM saw all other participants
		for (int i = 0; i < N + 1; i++)
		{
			System.out.println("=============== i = " + i + " with P = " + tsls[i].participants.size());
			Assert.assertEquals(tsls[i].participants.size(), N);
		}

		// disconnect
		for (int i = 0; i < N; i++)
		{
			peers[i].disconnect();
		}
	}
}
