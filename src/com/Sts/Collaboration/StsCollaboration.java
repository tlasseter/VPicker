package com.Sts.Collaboration;

import com.Sts.Actions.Wizards.Collaboration.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsCollaboration implements SessionListener
{
	StsModel model;
	StsCollaborationEntry collaborationEntry;
	StsCollaborationParticipant collaborationParticipant = null;
	boolean isLeader = false;

	// fields
	public boolean connected = false;
	public PeerChannel pc;
	public boolean created = false;
	public LeaderChannel lc;
	public boolean promoted = false;
	public boolean demoted = false;
	public Vector participants = new Vector();

	SessionManager sessionManager = null;
	Participant participant;

	static public StsCollaboration collaboration = null;

	static int leaderPort = 9001;
	static int nextPeerPort = 9011;
	static int lastPeerPort = nextPeerPort + 50;

	static boolean debug = false;

	private StsCollaboration(StsModel model, StsCollaborationEntry entry) throws Exception
	{
		this.model = model;
		this.collaborationEntry = entry;

		String username = System.getProperty("user.name");
		collaborationParticipant = collaborationEntry.getParticipantWithUsername(username);
		sessionManager = new SessionManager();
		sessionManager.addSessionListener(this);
		Runnable createSessionRunnable = new Runnable()
		{
			public void run()
			{
				isLeader = collaborationEntry.isLeader(collaborationParticipant);
				if (isLeader)
				{
					boolean isTest = collaborationEntry.name.equals("Test");
					if(!createSession(isTest))
					{
						if(!connectToSession(collaborationEntry, collaborationParticipant)) return;
					}
				}
				else
					if(!connectToSession(collaborationEntry, collaborationParticipant)) return;
			}
		};

		Thread createSessionThread = new Thread(createSessionRunnable);
		createSessionThread.start();
	}

	private boolean createSession(boolean isTest)
	{
		InetAddress address = null;
		try
		{
			address = InetAddress.getLocalHost();
			if (debug) System.out.println("Try to create session for leader: " + address.toString());
			sessionManager.createSession(address, leaderPort, "password");
			return true;
		}
		catch (Exception e) // couldn't connect on this port: try the next
		{
			if(!isTest)
			{
				StsMessageFiles.errorMessage("Failed to create session leader on " + address.toString() + ":" + leaderPort);
				StsException.outputException("StsCollaboration failed.", e, StsException.WARNING);
			}
			return false;
		}
	}

	private boolean connectToSession(StsCollaborationEntry collaborationEntry,
										StsCollaborationParticipant collaborationParticipant)
	{
		InetAddress leaderAddress = null;
		nextPeerPort = collaborationParticipant.minPort;
		lastPeerPort = collaborationParticipant.maxPort;
		int port = nextPeerPort;
		while (port <= lastPeerPort)
		{
			try
			{
				leaderAddress = InetAddress.getByName(collaborationEntry.getLeader().ipAddress);
				if (debug) System.out.println("Try to connect to session with leader: " + leaderAddress.toString());
				sessionManager.connectToSession(leaderAddress, leaderPort, port, "password");
				return true;
			}
			catch (Exception e) // couldn't connect on this port: try the next
			{
				if (port == lastPeerPort)
				{
					try
					{
						StsException.outputException("StsCollaboration failed. Peer couldn't connect to session.", e, StsException.WARNING);
						InetAddress address = InetAddress.getLocalHost();
						StsMessageFiles.errorMessage("Failed to create session peer on ports " + nextPeerPort + " - " + lastPeerPort + "\n" +
							"connecting to leader on " + leaderAddress.toString());
					}
					catch (Exception e1)
					{
						StsMessageFiles.errorMessage("Failed to create session peer on ports " + nextPeerPort + " - " + lastPeerPort + "\n" +
							"connecting to leader on " + leaderAddress.toString() + "\n" + "Could not define address for local host.");
					}
				}
                StsMessageFiles.errorMessage("Connection Error on port(" + port + "):" + e);
				port++;
			}
		}
		return false;
	}
/*
	private StsCollaboration(StsModel model) throws Exception
	{
		this.model = model;

		sessionManager = new SessionManager();
		sessionManager.addSessionListener(this);

		Runnable createSessionRunnable = new Runnable()
		{
			public void run()
			{
				try
				{
					if (mainDebug) System.out.println("Creating session for leader: " + InetAddress.getLocalHost().toString());
					sessionManager.createSession(InetAddress.getLocalHost(), leaderPort, "password");
				}
				catch (Exception e) // couldn't connect as leader, so try as peer
				{
					while (nextPeerPort <= lastPeerPort)
					{
						try
						{
							if (mainDebug) System.out.println("Connecting session to leader: " + InetAddress.getLocalHost().toString());
							sessionManager.connectToSession(InetAddress.getLocalHost(), leaderPort, nextPeerPort++, "password");
							break;
						}
						catch (Exception e1)
						{}
					}
				}
			}
		};

		Thread createSessionThread = new Thread(createSessionRunnable);
		createSessionThread.start();
	}
*/
	static public StsCollaboration construct(StsModel model, StsCollaborationEntry entry)
	{
		try
		{
			if(collaboration != null && collaboration.model == model && collaboration.collaborationEntry.equals(entry)) return collaboration;
			collaboration = new StsCollaboration(model, entry);
			return collaboration;
		}
		catch (Exception e)
		{
			return null;
		}
	}
/*
	static public StsCollaboration getInstance(StsModel model)
	{
		try
		{
			if (collaboration == null) collaboration = new StsCollaboration(model, null);
		}
		catch (Exception e)
		{
			return null;
		}
		return collaboration;
	}
*/
	static public StsCollaboration getCollaboration()
	{
        return collaboration;
	}

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
		try
		{
			if (pc == null)
			{
				StsException.systemError("StsCollaboration.constructor failed. No peerChannel created.");
				return;
			}
			participant = sessionManager.getIdentity();
			PeerObjectMsgListener ml = new PeerObjectMsgListener();
			pc.addListener(ml);
			System.out.println("Participant " + participant + " connected to session.");
			getDBFromLeader();
		}
		catch (Exception ex)
		{
			StsException.outputException("StsCollaboration.constructor() failed. Could create leader or peer connection.", ex, StsException.WARNING);
		}
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
		if (lc == null)
		{
			StsException.systemError("StsCollaboration.constructor failed. No leaderChannel created.");
			return;
		}
		participant = sessionManager.getIdentity();
		LeaderObjectMsgListener leaderMsgHandler = new LeaderObjectMsgListener(lc);
		lc.addListener(leaderMsgHandler); // setup inputStream to receive data from leader
		System.out.println("Leader " + participant + " created session.");
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
		PeerObjectMsgListener ml = new PeerObjectMsgListener();
		pc.addListener(ml);
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

	public void demoted()
	{
		try
		{
			Participant other = (Participant) participants.get(0);
			lc.demote(other);
		}
		catch (Exception e)
		{
			StsException.outputException("StsCollaboration.demoted() failed.", e, StsException.WARNING);
		}
	}

	public void promoted()
	{
		sessionManager.peerPromoted(participant);
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
		LeaderObjectMsgListener leaderMsgHandler = new LeaderObjectMsgListener(lc);
		lc.addListener(leaderMsgHandler); // setup inputStream to receive data from leader
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

	private boolean getDBFromLeader()
	{
		try
		{
			pc.sendToLeader(new RequestDBMessage(participant));
			model.getDatabase().close();
//			String pathname = System.getProperty("user.dirNo") + File.separator + "S2SCache" + File.separator + "db.temp";
//			StsFile file = StsFile.constructor(pathname);
//			if(!StsDBFile.fileOK(file)) return false;
			StsToolkit.runLaterOnEventThread(new Runnable()
			{
				public void run()
				{
					model.close();
					InputStream inputStream = pc.getInputStream();
					DataInputStream dis = new DataInputStream(inputStream);
					StsFile file = StsDBFile.copyToFile(dis);
					if (file == null)return;
					model = StsModel.constructor(file);
				}
			});

			return model != null;
		}
		catch (Exception e)
		{
			StsException.outputException("StsCollaboration.getDBFromLeader() failed.", e, StsException.WARNING);
			return false;
		}
	}

	public boolean broadcastTransaction(byte[] bytes, String description, byte transactionType)
	{
		try
		{
			if (debug && participants.size() > 0) System.out.println("Leader is sending " + bytes.length + " bytes to peers for transaction: " + description);
			TransactionMessage msg = new TransactionMessage(bytes, description, transactionType);
            if(debug) System.out.println("Leader:Message constructed.");
			lc.broadcast(null, msg);
            if(debug) System.out.println("Leader:Message broadcast.");
			return true;
		}
		catch (Exception e)
		{
			StsException.systemError("StsCollaboration.broadcastTransaction() failed.");
			return false;
		}
	}

	public boolean isLeader()
	{return sessionManager.isLeader();
	}

	public boolean hasPeers()
	{return sessionManager.hasPeers();
	}

	private class PeerObjectMsgListener implements PeerMessageListener
	{
		public LinkedList queue = new LinkedList();

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.MessageListener#messageReceived(java.io.Serializable)
		 */
		public void messageReceived(Serializable obj)
		{
			if (obj instanceof TransactionMessage)
			{
				TransactionMessage transactionMessage = (TransactionMessage) obj;
				byte[] transactionBytes = transactionMessage.transactionBytes;
				String description = transactionMessage.description;
				byte transactionType = transactionMessage.transactionType;
				if(debug) System.out.println("Peer " + participant + " received " + transactionBytes.length + " bytes from leader. Write transaction: " +
								   StsDBFile.transactionTypeString(transactionType));
				model.getDatabase().readTransaction(transactionBytes, description, transactionType, model);
			}
			else
				queue.addLast(obj);
		}
	}

	private class LeaderObjectMsgListener implements LeaderMessageListener
	{
		LeaderChannel lc;
		LeaderObjectMsgListener(LeaderChannel lc_)
		{
			lc = lc_;
		}

		public LinkedList queue = new LinkedList();

		/*
		 * (non-Javadoc)
		 *
		 * @see com.s2s.sm.MessageListener#messageReceived(java.io.Serializable)
		 */
		public void objectReceived(Participant p, Object obj)
		{
			if (obj instanceof RequestDBMessage)
			{
				try
				{
					OutputStream outputStream = lc.getOutputStream(p);
					model.getDatabase().writeToOutputStream(outputStream);
					//				lc.endBroadcast(p);
				}
				catch (Exception e)
				{
					System.out.println("TestLeaderMsgListener failed with error " + e.getMessage());
					e.printStackTrace();
				}
			}
			else
				queue.addLast(obj);
		}
	}

	static public void main(String[] args)
	{
		Properties props = System.getProperties();
		Enumeration names = props.propertyNames();
		while (names.hasMoreElements())
		{
			String name = (String) names.nextElement();
			String value = System.getProperty(name);
			System.out.println(name + " " + value);
		}
	}
}
