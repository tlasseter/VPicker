package com.Sts.Collaboration;

import com.Sts.Collaboration.protocol.*;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * This is the heart of the SessionManager (hence its name). The SessionManager
 * is used to connect to/create a session and to get a channel which is used to
 * interact with other peers in the system.
 *
 * @author druths
 */

public class SessionManager
{

	static final boolean debug = false;
	/*
	 * static SessionManager sessionManager = null;
	 *
	 * static public SessionManager getInstance() { if(sessionManager == null)
	 * sessionManager = new SessionManager(); return sessionManager; }
	 */
	// inner classes
	/**
	 * This class ensures that the local system is always connected to everybody
	 * in the session. When given a list of participants to connect to, it does so
	 * as quickly as possible.
	 *
	 * @author druths
	 */
	private class PeerListConnector extends Thread
	{

		// fields
		private LinkedList _unconnected_peers = new LinkedList();

		// constructors
		public PeerListConnector()
		{
			start();
		}

		// methods
		public synchronized void addUnconnectedPeers(Collection c)
		{
			if(debug) System.out.println("Adding " + c.size() + " peers");

			_unconnected_peers.addAll(c);

			notify();
		}

		public synchronized void addUnconnectedPeer(Participant p)
		{
			_unconnected_peers.addLast(p);

			notify();
		}

		public void run()
		{

			while (true)
			{

				synchronized (this)
				{
					if (_unconnected_peers.size() == 0)
					{
						try
						{
							if(debug) System.out.println("WAITING");
							wait();
						}
						catch (InterruptedException ie)
						{
							ie.printStackTrace();
							return;
						}

						if(debug) System.out.println("DONE WAITING");
					}
				}

				// get the next participant to connect to
				Participant p = (Participant) _unconnected_peers.removeFirst();

				// connect to the participant if necessary
				if (!_smanagers.containsKey(p) && !p.equals(_local_participant))
				{
					try
					{
						createConnection(p);
					}
					catch (IOException ioe)
					{
						ioe.printStackTrace();
						System.err.println("Unable to connect to peer " + p.getAddress() + ":" + p.getPort());
					}
				}
			}
		}
	}

	/**
	 * This class handles the loss of connections.
	 *
	 * @author druths
	 */
	private class SMListener implements SocketManagerListener
	{

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SocketManagerListener#socketManagerDeactivated(com.s2s.sm.SocketManager)
		 */
		public void socketManagerDeactivated(SocketManager sm)
		{


			Participant p = (Participant) _participants.get(sm);
			if(debug) System.out.println("SocketManagerDeactivated for " + p);

			if (p != null)
			{
				_smanagers.remove(p);
			}

			_participants.remove(sm);

			notifySessionListenersOfLostPeer(p);
		}

	}

	/**
	 * A container for session information.
	 *
	 * @author druths
	 */
	private class Session
	{

		private String _password;

		public void setPassword(String password)
		{
			_password = password;
		}

		public String getPassword()
		{
			return _password;
		}
	}

	/**
	 * This class handles incoming connections from other systems.
	 *
	 * @author druths
	 */
	private class SMSocketListener implements SessionServerSocketListener
	{

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionServerSocketListener#connectionFormed(com.s2s.sm.SessionServerSocket,
		 *      java.net.Socket)
		 */
		public void connectionFormed(SessionServerSocket sss, final Socket s)
		{

			synchronized (_promoting_lock)
			{
				if (_promoting)
				{
					System.err.println("Unable to accept participant " + s.getInetAddress() + ":" + s.getPort() + " due to promotion process");
					try
					{
						s.close();
					}
					catch (IOException ioe)
					{
						// we're not interested at ALL in what happens here.
					}

					return;
				}
			}

			// create a handler
			try
			{
				final SocketManager sm = new SocketManager(s);
                System.out.println("New Socket:" + s.toString());
				sm.addListener(new SMListener());

				sm.setObjectHandler(new SocketObjectHandler()
				{

					public void handleObject(SocketManager sm, Object obj)
					{
						Participant p = (Participant) obj;

						_smanagers.put(p, sm);
						_participants.put(sm, p);

						if (_is_leader == true)
						{
							_lc.addParticipant(p, sm);
						}

						if(debug) System.out.println(_local_participant + " added " + p + " to peer list");

						// send the set of systems the machine should connect to
						try
						{
							PeerList pl = new PeerList();
							pl.peers = new Vector(_smanagers.keySet());
							sm.writeObject(pl);
							if(debug) System.out.println("Sent peer list");
						}
						catch (IOException ioe)
						{
							System.err.println("Unable to send peer list to " + s.getInetAddress());

							sm.close();
						}

						notifySessionListenersOfNewPeer(p);

						return;
					}
				});
			}
			catch (IOException ioe)
			{
				System.err.println("Unable to create SocketManager for socket " + s.getInetAddress());
				return;
			}

		}

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.SessionServerSocketListener#sessionServerFailed(com.s2s.sm.SessionServerSocket)
		 */
		public void sessionServerFailed(SessionServerSocket sss)
		{
			// TODO Auto-generated method stub
		}
	}

	// fields
	private HashSet _s_listeners = new HashSet();

	private SessionServerSocket _sss;

	private Session _session;

	private Participant _local_participant;

	private boolean _is_leader = false;

	private LeaderChannel _lc;

	private PeerChannel _pc;

	private PeerListConnector _plm = new PeerListConnector();

	private Hashtable _smanagers = new Hashtable();

	private Hashtable _participants = new Hashtable();

	private Object _promoting_lock = new Object();

	private boolean _promoting = false;

	///// Methods
	public boolean isLeader()
	{return _is_leader;
	}

	public boolean hasPeers()
	{return _is_leader && _participants.size() > 0;
	}

	// listener methods
	/**
	 * Add an object that will be notified of session events.
	 */
	public void addSessionListener(SessionListener sl)
	{
		_s_listeners.add(sl);
	}

	/**
	 * Remove an object from the list of objects that are notified of session
	 * events.
	 */
	public void removeSessionListener(SessionListener sl)
	{
		_s_listeners.remove(sl);
	}

	private void notifySessionListenersOfCreated(LeaderChannel lc)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.createdSession(this, lc);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".connectedToSession(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfConnected(PeerChannel pc)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.connectedToSession(this, pc);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".connectedToSession(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfNewPeer(Participant p)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.connectedToNewPeer(this, p);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".connectedToNewPeer(" + this +", " + p + ") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfLostPeer(Participant p)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.lostConnectionToPeer(this, p);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".lostConnectionToPeer(" + this +", " + p + ") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfDisconnected()
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.disconnectedFromSession(this);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".disconnectedFromSession(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfDemotion(PeerChannel pc)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.demoted(this, pc);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".demoted(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfPromotion(LeaderChannel lc)
	{
		if(_s_listeners.size() == 0) System.err.println("No sessionListers exist to notify of promotion.");
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.promoted(this, lc);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".promoted(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	private void notifySessionListenersOfPeerPromoted(Participant p)
	{
		Iterator i = _s_listeners.iterator();
		SessionListener sl;
		while (i.hasNext())
		{
			sl = (SessionListener) i.next();
			try
			{
				sl.peerPromoted(this, p);
			}
			catch (RuntimeException rte)
			{
				System.err.println("ERROR: " + sl + ".peerPromoted(" + this +") threw exception");
				rte.printStackTrace();
			}
		}
	}

	// connection methods
	/**
	 * This method creates a local session including only this machine.
	 *
	 * @param addr
	 *          is the local address that other peers can use to join this
	 *          session.
	 * @param port
	 *          is the local port that other peers will connect on.
	 * @param password
	 *          is the password that other peers must provide to be admitted.
	 */
	public synchronized void createSession(InetAddress addr, int port, String password) throws AlreadyInSessionException, IOException
	{

		if (_session != null)
		{
			throw new AlreadyInSessionException();
		}

		_local_participant = new Participant(addr, port);

		// create a session server socket
		_sss = new SessionServerSocket(addr, port);
		_sss.addListener(new SMSocketListener());

		// create a leader channel
		_lc = new LeaderChannel(this, _smanagers);
		_is_leader = true;

		_session = new Session();
		_session.setPassword(password);

		// Notify listeners
		notifySessionListenersOfCreated(_lc);
	}

	/**
	 * This method makes this SessionManager join an already existing session at
	 * the remote address/port specified. The password is used to authenticate
	 * with the remote system. The remote system must be a host in the
	 * collaboration session.
	 */
	public synchronized void connectToSession(InetAddress addr, int port, int local_port, final String password) throws AlreadyInSessionException, IOException
	{

		if (_session != null)
		{
			throw new AlreadyInSessionException();
		}

		_local_participant = new Participant(InetAddress.getLocalHost(), local_port);

		SocketManager sm = createConnection(addr, port);

		// create a session server socket
		_sss = new SessionServerSocket(InetAddress.getLocalHost(), local_port);
		_sss.addListener(new SMSocketListener());

		// get the machines to connect to
		sm.setObjectHandler(new SocketObjectHandler()
		{
			public void handleObject(SocketManager sm, Object obj)
			{

//				System.out.println(this +" got obj = " + obj);

				PeerList peerList = null;
				peerList = (PeerList) obj;

				if(debug) System.out.println("Got peer list");

				// start the process of connecting to the other peers
				_plm.addUnconnectedPeers(peerList.peers);

				// create a peer channel
				_pc = new PeerChannel(SessionManager.this, sm);

				_session = new Session();
				_session.setPassword(password);

				// Do any necessary initialization
				// see my notes in createSession
				// create a connection to the leader (addr, port)
				notifySessionListenersOfConnected(_pc);
			}
		});
	}

	private SocketManager createConnection(Participant p) throws IOException
	{

		SocketManager sm = new SocketManager(p.getAddress(), p.getPort());
		sm.addListener(new SMListener());

		// send the new system our identity.
		sm.writeObject(_local_participant);

		// put info about this new participant in an easy to find place.
		_smanagers.put(p, sm);
		_participants.put(sm, p);

		System.out.println("Connection created to " + p);

		notifySessionListenersOfNewPeer(p);

		return sm;
	}

	private SocketManager createConnection(InetAddress addr, int port) throws IOException
	{

		Participant p = new Participant(addr, port);

		return createConnection(p);
	}

	/**
	 * Return the complete list of currently connected remote systems and the
	 * local system. If the session manager is not in a session, then the set is
	 * empty.
	 */
	public Set getParticipants()
	{
		return _smanagers.keySet();
	}

	public synchronized void disconnect() throws StillLeaderException
	{

		if (_session == null)
		{

			return;
		}

		// close the channel
		if (_is_leader)
		{

			if (_smanagers.size() > 0)
			{

				throw new StillLeaderException();
			}

			_lc.deactivate();
		}
		else
		{

			_pc.disconnect();
		}

		// close the server socket if it exists
		if (_sss != null)
		{

			_sss.close();
		}

		// close all active socket managers
		Iterator i = new Vector(_smanagers.keySet()).iterator();

		while (i.hasNext())
		{

			 ( (SocketManager) _smanagers.get(i.next())).close();
		}

		// done!
	}

	public Participant getIdentity()
	{
		return _local_participant;
	}

	void setPromotionLock(boolean promoting)
	{

		synchronized (_promoting_lock)
		{
			_promoting = promoting;
		}
	}

	void peerPromoted(Participant p)
	{

		System.out.println("PEER PROMOTED");

		if (p.equals(_local_participant))
		{
			System.out.println("PROMOTING!!");

			_is_leader = true;
			_lc = new LeaderChannel(this, _smanagers);

			notifySessionListenersOfPromotion(_lc);

			// make the former PeerChannel unusable
			_pc.deactivate();
		}
		else
		{
			_pc.swapLeaders( (SocketManager) _smanagers.get(p));
			notifySessionListenersOfPeerPromoted(p);
		}
	}

	void demote(Participant p)
	{

		_is_leader = false;
		_pc = new PeerChannel(this, (SocketManager) _smanagers.get(p));

		notifySessionListenersOfDemotion(_pc);
	}
}
