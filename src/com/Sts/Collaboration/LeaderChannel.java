package com.Sts.Collaboration;

import com.Sts.Collaboration.protocol.*;

import java.io.*;
import java.util.*;

/**
 * This class defines the interface used to interact with peers when the local
 * system has leader status.
 *
 * @author druths
 */
public class LeaderChannel implements Channel
{

	// constants
	public static final int PROMOTION_TIMEOUT = 10000; // 10 seconds

	// inner classes
	private class LeaderObjectHandler implements SocketObjectHandler
	{

		private Participant _p;

		public LeaderObjectHandler(Participant p)
		{
			_p = p;
		}

		public void handleObject(SocketManager sm, Object obj)
		{

			if (obj instanceof PeerList)
			{
				if (_in_promotion_process)
				{
					completePromotionProcess(_p, (PeerList) obj);
				}
			}

			notifyListenersOfObjectReceived(_p, obj);
		}
	}

	private class ByteBundler extends OutputStream
	{

		/*
		 * if singleReceiver not null, broadcast to all participants; otherwise only
		 * to this singleReceiver
		 */
		Participant singleReceiver = null;

		/*
		 * (non-Javadoc)
		 *
		 * @see java.io.OutputStream#close()
		 */public void close() throws IOException
		{
			throw new IOException("The OuputStream cannot be closed on the LeaderChannel");
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.io.OutputStream#flush()
		 */
		public void flush() throws IOException
		{
			return;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.io.OutputStream#write(byte[], int, int)
		 */
		public void write(byte[] arg0, int arg1, int arg2) throws IOException
		{
			if (_active == false)
			{
				throw new IOException("Leader not active");
			}
			byte[] array = new byte[arg2 - arg1];
			for (int i = 0; i < (arg2 - arg1); i++)
			{
				array[i] = arg0[arg1 + i];
			}
			broadcast(singleReceiver, new ByteBundle(array));
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.io.OutputStream#write(byte[])
		 */
		public void write(byte[] arg0) throws IOException
		{
			if (_active == false)
			{
				throw new IOException("Leader not active");
			}
			broadcast(singleReceiver, new ByteBundle(arg0));
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.io.OutputStream#write(int)
		 */
		public void write(int arg0) throws IOException
		{
			if (_active == false)
			{
				throw new IOException("Leader not active");
			}
			broadcast(singleReceiver, new ByteBundle(arg0));
		}
	}

	/**
	 * This object enforces a timelimit on the promotion process.
	 *
	 * @author druths
	 */
	private class PromotionTimerTask extends TimerTask
	{

		/*
		 * (non-Javadoc)
		 *
		 * @see java.util.TimerTask#run()
		 */
		public void run()
		{
			failPromotionProcess();
		}
	}

	// fields
	private Hashtable _peers;

	private SessionManager _sm;

	private boolean _active = true;

	private Object _promotion_lock = new Object();

	private boolean _in_promotion_process = true;

	private ByteBundler _byte_bundler = new ByteBundler();

	private HashSet _listeners = new HashSet();

	// constructors
	public LeaderChannel(SessionManager sm, Hashtable peers)
	{
		_sm = sm;
		_peers = peers;
		addParticipants(peers);
	}

	//// methods
	private void addParticipants(Hashtable peers)
	{
		Iterator i = peers.keySet().iterator();
		while (i.hasNext())
		{
			Participant p = (Participant) i.next();
			SocketManager sm = (SocketManager) peers.get(p);
			addParticipant(p, sm);
		}
	}

	public void addListener(LeaderMessageListener lcl)
	{
		_listeners.add(lcl);
	}

	public void removeListener(LeaderMessageListener lcl)
	{
		_listeners.remove(lcl);
	}

	private void notifyListenersOfObjectReceived(Participant p, Object obj)
	{
		Iterator i = _listeners.iterator();
		while (i.hasNext())
		{
			 ( (LeaderMessageListener) i.next()).objectReceived(p, obj);
		}
		return;
	}

	public void addParticipant(Participant p, SocketManager sm)
	{
		sm.setObjectHandler(new LeaderObjectHandler(p));
	}

	/**
	 * Send a message to a specific peer.
	 *
	 * @throws IOException
	 *           if the message can't be sent.
	 */
	public void sendToPeer(Participant p, Serializable msg) throws IOException
	{
		if (!_active)
		{
			throw new IOException("Deactivated LeaderChannel");
		}
		inner_sendToPeer(p, msg);
	}

	private void inner_sendToPeer(Participant p, Serializable msg) throws IOException
	{
		SocketManager sm = (SocketManager) _peers.get(p);
		if (sm == null)
		{
			throw new IOException("Peer " + p + " was not found");
		}
		if (SessionManager.debug)
		{
			System.out.println("inner_sendToPeer " + p + " message " + msg.getClass().getName());
		}
		sm.writeObject(msg);
	}

	private void inner_endBroadcast(Participant p) throws IOException
	{
		SocketManager sm = (SocketManager) _peers.get(p);
		if (sm == null)
		{
			throw new IOException("Peer " + p + " was not found");
		}
		sm.shutdownInput();
	}

	/**
	 * Sends a message to all participants unless singleReceiver is not null in
	 * which case send to just this receiver. The IOException is thrown late -
	 * meaning that it isn't thrown until the message is sent to as many
	 * participants as possible. The IOException is of type that will contain all
	 * the errors that occured.
	 */
	public void broadcast(Participant singleReceiver, Serializable msg) throws IOException
	{
		if (!_active)
		{
			throw new IOException("Deactivated LeaderChannel");
		}
		inner_broadcast(singleReceiver, msg);
	}

	void inner_broadcast(Participant singleReceiver, Serializable msg) throws IOException
	{
		if (singleReceiver != null)
		{
			try
			{
				if(SessionManager.debug)
				{
					System.out.println("Sending to " + singleReceiver);
				}
				inner_sendToPeer(singleReceiver, msg);
			}
			catch (IOException ioe)
			{
				System.err.println("inner_broadcast failed for single receiver: " + singleReceiver.toString());
				System.err.println("\t" + ioe.getMessage());
			}
		}
		else
		{
			IOException saved = null;

			Iterator i = _peers.keySet().iterator();
			while (i.hasNext())
			{
				Participant p = (Participant) i.next();
				try
				{
					if (SessionManager.debug)
					{
						System.out.println("Sending to " + p);
					}
					inner_sendToPeer(p, msg);
				}
				catch (IOException ioe)
				{
					saved = ioe;
				}
			}
			if (saved != null)
			{
				System.err.println("Unable to reach all peers in list");
				System.err.println("\t" + saved.getMessage());
			}
		}
	}

	public void endBroadcast(Participant singleReceiver)
	{
		if (singleReceiver != null)
		{
			try
			{
				if (SessionManager.debug)
				{
					System.out.println("Sending to " + singleReceiver);
				}
				inner_endBroadcast(singleReceiver);
			}
			catch (Exception e)
			{
				System.err.println("endBroadcast failed for single receiver: " + singleReceiver.toString());
				System.err.println("\t" + e.getMessage());
			}
		}
		else
		{
			Exception saved = null;

			Iterator i = _peers.keySet().iterator();
			while (i.hasNext())
			{
				Participant p = (Participant) i.next();
				try
				{
					if (SessionManager.debug)
					{
						System.out.println("Sending to " + p);
					}
					inner_endBroadcast(p);
				}
				catch (Exception e)
				{
					saved = e;
				}
			}
			if (saved != null)
			{
				System.err.println("Unable to end broadcast to all peers in list");
				System.err.println("\t" + saved.getMessage());
			}
		}
	}

	/**
	 * A call to this method causes the local system to attempt to give the
	 * specified participant the leader status. The success of this procedure is
	 * not guaranteed and is non-blocking. SessionListeners will be notified of
	 * success. During the time a promotion event is occurring, the leader does
	 * not accept any new participants and messages cannot be sent.
	 *
	 * @param new_leader
	 * @throws IOException
	 */
	public void demote(Participant new_leader) throws IOException
	{
		if (!_active)
		{
			throw new IOException("Deactivated LeaderChannel");
		}
		synchronized (_promotion_lock)
		{
			deactivate();
			// keep activities that shouldn't happen during a promotion process
			// from happening
			_in_promotion_process = true;
			_sm.setPromotionLock(true);
			// set a timer to cause an automatic fail if the promotion process
			// takes too long
			new Timer().schedule(new PromotionTimerTask(), PROMOTION_TIMEOUT);
			try
			{
				inner_sendToPeer(new_leader, new RequestPeerListMsg());
				System.out.println("Demoted " + _sm.getIdentity() + "  Promoting " + new_leader);
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
				failPromotionProcess();
				return;
			}
			//new PeerCommReader(new_leader, (SocketManager)_peers.get(new_leader));
		}
	}

	/**
	 * Once the peer has accepted the new leader status, the current leader
	 * notified everybody that it is the new leader.
	 */
	private void completePromotionProcess(Participant p, PeerList pl)
	{
		// check if the peer lists match
		HashSet pl1 = new HashSet(pl.peers);
		pl1.add(p);
		HashSet pl2 = new HashSet(_peers.keySet());
		pl2.add(_sm.getIdentity());
		Iterator i = pl1.iterator();
		while (i.hasNext())
		{
			pl2.remove(i.next());
		}
		if (pl2.size() > 0)
		{
			System.err.println("Peer Lists are not the same");
			failPromotionProcess();
		}
		else
		{
			if(SessionManager.debug) System.out.println("Peer lists are the same.  Broadcasting");
			try
			{
				inner_broadcast(null, new PromoteMsg(p));
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
			}
			_sm.demote(p);
		}
	}

	/**
	 * This method is called to cancel the promotion process should anything go
	 * wrong during the process.
	 */
	private void failPromotionProcess()
	{
		synchronized (_promotion_lock)
		{
			if (_in_promotion_process)
			{
				activate();
				// reset promotion block
				_sm.setPromotionLock(false);
				_in_promotion_process = false;
			}
		}
		// TODO add notification here?
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see com.s2s.sm.Channel#isActive()
	 */
	public boolean isActive()
	{
		return _active;
	}

	/**
	 * This method temporarily makes this leader channel unusable.
	 *
	 */
	void deactivate()
	{
		_active = false;
	}

	void activate()
	{
		_active = true;
	}

	public OutputStream getOutputStream()
	{
		_byte_bundler.singleReceiver = null;
		return _byte_bundler;
	}

	public OutputStream getOutputStream(Participant singleReceiver)
	{
		_byte_bundler.singleReceiver = singleReceiver;
		return _byte_bundler;
	}

}
