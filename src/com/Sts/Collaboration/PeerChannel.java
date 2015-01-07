package com.Sts.Collaboration;

import com.Sts.Collaboration.protocol.*;

import java.io.*;
import java.util.*;

/**
 * This class provides the ability to receive messages from the leader.
 */
public class PeerChannel implements Channel
{

	// inner classes
	/**
	 * This class is used during a possible promotion event. It spawns a new
	 * thread to send the PeerList back to the leader.
	 */
	private class PeerListSender extends Thread
	{

		public void run()
		{
			try
			{
				PeerList pl = new PeerList();
				pl.peers = new Vector(_sm.getParticipants());
				_leader.writeObject(pl);
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
				disconnect();
			}
		}
	}

	private class SendingThread extends Thread
	{

		private List _pending = new LinkedList();

		public void addPending(Object obj)
		{
			synchronized (this)
			{
				_pending.add(obj);
				this.notify();
			}
		}

		public void run()
		{
			Object obj;
			try
			{
				while (true)
				{
					synchronized (this)
					{
						if (_pending.isEmpty())
						{
							this.wait();
						}
						obj = _pending.remove(0);
					}
					_leader.writeObject(obj);
				}
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
				disconnect();
			}
			catch (InterruptedException ie)
			{
				ie.printStackTrace();
			}
		}
	}

	/**
	 * This is a thread that performs retrieval and handling of messages sent by
	 * the leader.
	 */
	private class PeerObjectHandler implements SocketObjectHandler
	{
		public void handleObject(SocketManager sm, Object msg)
		{
			// handle the beginning of a promotion/demotion process
			if (msg instanceof RequestPeerListMsg)
			{
				// start a new thread to send the peer list
				if(SessionManager.debug) System.out.println("New leader " + _sm.getIdentity() + "received RequestPeerListMsg. Sending peer list.");
				new PeerListSender().start();
				return;
			}
			// handle a promotion event
			if (msg instanceof PromoteMsg)
			{
				if(SessionManager.debug) System.out.println("PromoteMsg received by " + _sm.getIdentity() + ". New leader is " + ((PromoteMsg) msg).new_leader);
				_sm.peerPromoted( ( (PromoteMsg) msg).new_leader);
				return;
			}
			if (msg instanceof ByteBundle)
			{
				ByteBundle bb = (ByteBundle) msg;
				try
				{
					if (bb.hasInt)
					{
						_pos.write(bb.val);
					}
					else
					{
						_pos.write(bb.array);
					}
				}
				catch (IOException ioe)
				{
					ioe.printStackTrace();
				}

				return;
			}
			// if it wasn't one of the above, pass handling off to the
			// MessageListener.
			synchronized (PeerChannel.this)
			{
				// queue up messages until we get a MessageListener to give these
				// to.
				if (_ml == null)
				{
					_msg_queue.addLast(msg);
				}
				else
				{
					_ml.messageReceived( (Serializable) msg);
				}
			}
		}
	}

	// fields
	private boolean _active = true;

	private SocketManager _leader;

	private PeerObjectHandler _poh;

	private SendingThread _sender;

	private PeerMessageListener _ml = null;

	private LinkedList _msg_queue = new LinkedList();

	private SessionManager _sm;

	private PipedInputStream _pis;

	private PipedOutputStream _pos;

	// constructors
	public PeerChannel(SessionManager sm, SocketManager leader)
	{
		_leader = leader;
		_sm = sm;
		_sender = new SendingThread();
		_sender.start();

		try
		{
			_pis = new PipedInputStream();
			_pos = new PipedOutputStream(_pis);
		}
		catch (IOException ioe)
		{
			ioe.printStackTrace();
			System.err.println("This shouldn't happen");
		}
		_poh = new PeerObjectHandler();
		_leader.setObjectHandler(_poh);
	}

	//// methods
	/**
	 * This method changes out leaders that this channel gets messages from.
	 */
	public void swapLeaders(SocketManager leader)
	{
		// replace the handler
		_leader.setObjectHandler(null);
		_leader = leader;
		_leader.setObjectHandler(_poh);

		// stop the last one
		_sender.interrupt();

		// start a new one
		_sender = new SendingThread();
		_sender.start();
	}

	/**
	 * Disconnect the local system from the session.
	 */
	public void disconnect()
	{
		_active = false;
		_leader.setObjectHandler(null);
		_leader.close();
		try
		{
			_pos.close();
		}
		catch (IOException ioe)
		{
			// nothing to do here, doesn't matter
		}
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
	 * This method deactivates the local channel without closing the underlying
	 * connection to the leader.
	 */
	public void deactivate()
	{
		_active = false;
		_leader.setObjectHandler(null);
	}

	/**
	 * Set the object that will receive notifications when messages arrive. If
	 * messages have already arrived, then the object will be notified of each in
	 * the order they were received.
	 *
	 * @param ml
	 *          the object to be notified.
	 */
	public synchronized void addListener(PeerMessageListener ml)
	{
		_ml = ml;
		dequeueMessages();
	}

	/**
	 * Notify the MessageListener of all the messages that were received before it
	 * was registered.
	 */
	private void dequeueMessages()
	{
		while (_msg_queue.size() > 0)
		{
			try
			{
				_ml.messageReceived( (Serializable) _msg_queue.removeFirst());
			}
			catch (RuntimeException rte)
			{
				rte.printStackTrace();
			}
		}
	}

	public InputStream getInputStream()
	{
		return _pis;
	}

	public void sendToLeader(Serializable s) throws IOException
	{
		_sender.addPending(s);
	}
}
