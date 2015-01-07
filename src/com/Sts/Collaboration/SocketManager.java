package com.Sts.Collaboration;

import java.io.*;
import java.net.*;
import java.util.*;

/**
 * The socket manager manages the liveness of a given socket as well as the
 * serialization and deserialization of objects being sent through the streams.
 *
 * @author druths
 */
class SocketManager
{

	// fields
	private Socket _s;

	private ObjectOutputStream _oos;

	private ObjectInputStream _ois;

	private HashSet _listeners = new HashSet();

	private Thread _reader;

	private SocketObjectHandler _soh = null;

	public SocketManager(InetAddress addr, int port) throws IOException
	{
		_s = new Socket(addr, port);
		_oos = new ObjectOutputStream(_s.getOutputStream());
		startReader();
	}

	public SocketManager(Socket s) throws IOException
	{
		_s = s;
		_oos = new ObjectOutputStream(_s.getOutputStream());
		startReader();
	}

	private void startReader()
	{
		Thread t = new Thread()
		{

			public void run()
			{

				Object obj;

				while (true)
				{
					try
					{
						obj = readObject();

						synchronized (SocketManager.this)
						{
							if (_soh == null)
							{
								SocketManager.this.wait();
							}
						}
						_soh.handleObject(SocketManager.this, obj);
					}
					catch (IOException ioe)
					{
//						ioe.printStackTrace();
						return;
					}
					catch (ClassNotFoundException cnfe)
					{
						cnfe.printStackTrace();
						return;
					}
					catch (InterruptedException ie)
					{
						ie.printStackTrace();
						return;
					}
				}
			}
		};

		t.start();
		_reader = t;
	}

	public synchronized void setObjectHandler(SocketObjectHandler soh)
	{
		_soh = soh;
		notifyAll();
	}

	private Object readObject() throws IOException, ClassNotFoundException
	{
		if (_ois == null)
		{
            System.out.println("Creating reader..." + _s.toString());
			_ois = new ObjectInputStream(_s.getInputStream());
		}

		try
		{
			Object obj = _ois.readObject();
            return obj;
		}
		catch (IOException ioe)
		{
			close();
			throw ioe;
		}
		catch (ClassNotFoundException cnfe)
		{
			close();
			throw cnfe;
		}
	}

	public void writeObject(Object obj) throws IOException
	{
		try
		{

			_oos.writeObject(obj);
            _oos.flush();
		}
		catch (IOException ioe)
		{
			close();
			throw ioe;
		}
	}

	public void shutdownInput() throws IOException
	{
		_s.shutdownInput();
	}

	public void close()
	{
		try
		{
			_s.close();
		}
		catch (IOException ioe)
		{
			// ok... the socket may already be closed.
		}

		notifyListenersOfDeactivation();
	}

	public void addListener(SocketManagerListener sml)
	{

		_listeners.add(sml);
	}

	public void removeListener(SocketManagerListener sml)
	{

		_listeners.remove(sml);
	}

	private void notifyListenersOfDeactivation()
	{

		Iterator i = _listeners.iterator();

		while (i.hasNext())
		{

			SocketManagerListener sml = (SocketManagerListener) i.next();

			try
			{

				sml.socketManagerDeactivated(this);
			}
			catch (RuntimeException rte)
			{

				rte.printStackTrace();
			}
		}
	}
}
