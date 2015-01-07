package com.Sts.IO;

import java.io.*;
import java.util.*;
/**
 This reader/writer lock prevents reads from occurring while writes are
 in progress, and it also prevents multiple writes from
 happening simultaneously. Multiple read operations can run in parallel,
 however. Reads take priority over writes, so any read operations that
 are pending while a write is in progress will execute before any
 subsequent writes execute. Writes are guaranteed to execute in the
 order in which they were requested -- the oldest request is processed
 first.

 You should use the lock as follows:

 public class Data_structure_or_resource
 {
  Reader_writer lock = new Reader_writer();

  public void access( )
  {
   try
   {   lock.request_read();
  // do the read/access operation here.
   }
   finally
   {   lock.read_accomplished();
   }
  }

  public void modify( )
  {
   try
   {   lock.request_write();
  // do the write/modify operation here.
   }
   finally
   {   lock.write_accomplished();
   }
  }
 }

 This implementation is based on the one in Doug Lea's Concurrent
 Programming in Java (Addison Wesley, 1997, pp. 300-303), I've
 simplified the code (and cleaned it up) and added the nonblocking
 acquisition methods.  I've also made the lock a standalone class rather
 than a base class from which you have to derive. You might also want to
 look at the very different implementation of the reader/writer lock in
 Scott Oaks and Henry Wong's Java Threads (O'Reilly, 1997, pp.
 180-187).

 @author Allen I. Holub

 */

public class StsIOLock
{
	private int activeReaders; // = 0
	private int waitingReaders; // = 0
	private int activeWriters; // = 0
	private final boolean debug = false;

	/******************************************************************
	 I keep a linked list of writers waiting for access so that I can
	 release them in the order that the requests were received.  The size of
	 this list is the "waiting writers" count.  Note that the monitor of the
	 Reader_writer object itself is used to lock out readers
	 while writes are in progress, thus there's no need for a separate
	 "reader_lock."

	 */

	private final LinkedList writerLocks = new LinkedList();

	/******************************************************************
	 Request the read lock. Block until a read operation can be performed
	 safely.  This call must be followed by a call to
	 read_accomplished() when the read operation completes.

	 */

	public synchronized void requestRead()
	{
		requestRead("");
	}

	public synchronized void requestRead(String requestor)
	{
		if (activeWriters == 0 && writerLocks.size() == 0)
		{
			if(debug) System.out.println(requestor + " requesting read: ok.");
			++activeReaders;
		}
		else
		{
			++waitingReaders;
			try
			{
				if(debug) System.out.println(requestor + " requesting read: not ok. Waiting...");
				wait();
			}
			catch (InterruptedException e)
			{}
		}
	}

	/******************************************************************
	 This version of read() requests read access and returns
	 true if you get it. If it returns false, you may not
	 safely read from the guarded resource. If it returns true, you
	 should do the read, then call read_accomplished in the
	 normal way. Here's an example:

	 public void read()
	 {   if( lock.request_immediate_read() )
	  {   try
	   {
	 // do the read operation here
	   }
	   finally
	   {   lock.read_accomplished();
	   }
	  }
	  else
	   // couldn't read safely.
	 }
	 */

	public synchronized boolean requestImmediateRead()
	{
		if (activeWriters == 0 && writerLocks.size() == 0)
		{
			++activeReaders;
			return true;
		}
		return false;
	}

	/******************************************************************      |
	 Release the lock. You must call this method when you're done
	 with the read operation.
	 */


	public synchronized void readAccomplished()
	{
		readAccomplished("");
	}

	public synchronized void readAccomplished(String notifier)
	{
		if (--activeReaders == 0)
		{
			if(debug) System.out.println(notifier + "notifying read accomplished.");
			notifyWriters(notifier);
		}
	}

	/******************************************************************      |
	 Request the write lock. Block until a write operation can be performed
	 safely. Write requests are guaranteed to be executed in the order
	 received. Pending read requests take precedence over all write
	 requests.  This call must be followed by a call to
	 write_accomplished() when the write operation completes.

	 */

	public synchronized void requestWrite()
	{
		requestWrite("");
	}

	public void requestWrite(String requestor)
	{
		// This method can't be synchronized or there'd be a nested-monitor
		// lockout problem: We have to acquire the lock for "this" in
		// order to modify the fields, but that lock must be released
		// before we start waiting for a safe time to do the writing.
		// If request_write() were synchronized, we'd be holding
		// the monitor on the Reader_writer lock object while we were
		// waiting. Since the only way to be released from the wait is
		// for someone to call either read_accomplished()
		// or write_accomplished() (both of which are synchronized),
		// there would be no way for the wait to terminate.

		Object lock = new Object();
		synchronized (lock)
		{
			synchronized (this)
			{
				boolean okay_to_write = writerLocks.size() == 0 && activeReaders == 0 && activeWriters == 0;
				if (okay_to_write)
				{
					if(debug) System.out.println(requestor + " requesting write: ok.");
					++activeWriters;
					return; // the "return" jumps over the "wait" call
				}
				writerLocks.addLast(lock);
			}
			try
			{
				if(debug) System.out.println(requestor + " requesting write: not ok. Waiting...");
				lock.wait();
			}
			catch (InterruptedException e)
			{}
		}
	}

	/******************************************************************
	 This version of the write request returns false immediately
	 (without blocking) if any read or write operations are in progress and
	 a write isn't safe; otherwise, it returns true and acquires the
	 resource. Use it like this:


	 public void write()
	 {   if( lock.request_immediate_write() )
	  {   try
	   {
	 // do the write operation here
	   }
	   finally
	   {   lock.write_accomplished();
	   }
	  }
	  else
	   // couldn't write safely.
	 }
	 */

	synchronized public boolean requestImmediateWrite()
	{
		if (writerLocks.size() == 0 && activeReaders == 0 && activeWriters == 0)
		{++activeWriters; return true;
		}
		return false;
	}

	/******************************************************************
	 Release the lock. You must call this method when you're done
	 with the read operation.
	 */

	public synchronized void writeAccomplished()
	{
		writeAccomplished("");
	}

	public synchronized void writeAccomplished(String notifier)
	{
		// The logic here is more complicated than it appears.
		// If readers have priority, you'll  notify them. As they
		// finish up, they'll call read_accomplished(), one at
		// a time. When they're all done, read_accomplished() will
		// notify the next writer. If no readers are waiting, then
		// just notify the writer directly.

		--activeWriters;
		if(debug) System.out.println(notifier + " notifying write accomplished.");
		if (waitingReaders > 0) // priority to waiting readers
		{
			notifyReaders(notifier);
		}
		else
		{
			notifyWriters(notifier);
		}
	}

	/******************************************************************
	 Notify all the threads that have been waiting to read.
	 */

	public synchronized void notifyReaders()
	{
		notifyReaders("");
	}

	private void notifyReaders(String notifier) // must be accessed from a
	{ //  synchronized method
		activeReaders += waitingReaders;
		waitingReaders = 0;
		if(debug) System.out.println(notifier + " notifying waiting readers...");
		notifyAll();
	}

	/******************************************************************
	 Notify the writing thread that has been waiting the longest.
	 */

	public synchronized void notifyWriters()
	{
		notifyWriters("");
	}

	private void notifyWriters(String notifier) // must be accessed from a
	{ //  synchronized method
		if (writerLocks.size() > 0)
		{
			if(debug && activeWriters > 0) System.out.println(notifier + " notifying the oldest waiting writer...");
			Object oldest = writerLocks.removeFirst();
			++activeWriters;
			synchronized (oldest)
			{
				oldest.notify();
			}
		}
	}

	/*******************************************************************
	 The Test class is a unit test for the other code in the current file. Run the test with:

	 java com.holub.asynch.Reader_writer\$Test

	 (the backslash isn't required with windows boxes), and don't include
	 this class file in your final distribution.  The output could
	 vary in trivial ways, depending on system timing. The read/write order
	 should be exactly the same as in the following sample:

	 Starting w/0
	  w/0 writing
	 Starting r/1
	 Starting w/1
	 Starting w/2
	 Starting r/2
	 Starting r/3
	  w/0 done
	 Stopping w/0
	  r/1 reading
	  r/2 reading
	  r/3 reading
	  r/1 done
	 Stopping r/1
	  r/2 done
	  r/3 done
	 Stopping r/2
	 Stopping r/3
	  w/1 writing
	  w/1 done
	 Stopping w/1
	  w/2 writing
	  w/2 done
	 Stopping w/2
	 */

	static public void main(String[] args)
	{
		String filename = System.getProperty("user.dirNo") + File.separator + "testFile";
		Test t = new Test(filename);
	}

	public static class Test
	{
//		Resource resource = new Resource();
		/**
		 Test by creating several readers and writers. The initial write
		 operation (w/0) should complete before the first read (r/1) runs. Since
		 readers have priority, r/2 and r/3 should run before w/1; and r/1, r/2
		 and r3 should all run in parallel.  When all three reads complete, w1
		 and w2 should execute sequentially in that order.
		 */

		public Test(String filename)
		{
/*
			if (!resource.readIfPossible())
			{
				System.out.println("Immediate read request didn't work");
			}
			if (!resource.writeIfPossible())
			{
				System.out.println("Immediate write request didn't work");
			}
*/
		    File file = new File(filename);
			long length = 0;
			if(file.exists()) length = file.length();
			System.out.println("File " + filename + " has " + length + " bytes.");

	        StsIOLock lock = new StsIOLock();
		    new Writer(lock, "w/0", filename, 1000).start();
			new Reader(lock, "r/1", filename).start();
			new Writer(lock, "w/1", filename, 500).start();
			new Writer(lock, "w/2", filename, 250).start();
			new Reader(lock, "r/2", filename).start();
			new Reader(lock, "r/3", filename).start();
		}

		/** A reader thread. Reads a lock object, name descriptor, and filename */
		class Reader extends Thread
		{
			private String name;
			private String filename;
			FileInputStream fis;
			StsIOLock lock;

			Reader(StsIOLock lock, String name, String filename)
			{
				this.lock = lock;
				this.name = name;
				this.filename = filename;
				try
				{
					fis = new FileInputStream(filename);
				}
				catch(Exception e) { }
			}

			public void run()
			{
				System.out.println("Starting " + name);
				long nBytesTotal = 0;
				byte[] buf = new byte[1024];
				try
				{
					lock.requestRead("test read");
					while(true)
					{
						int nBytesRead = fis.read(buf);
						if (nBytesRead > 0) nBytesTotal += nBytesRead;
						else
						{
							lock.readAccomplished("test read");
							System.out.println("Stopping " + name + " read " + nBytesTotal + " bytes.");
							break;
						}
					}
					Thread.currentThread().sleep(100);
				}
				catch(Exception e) {}
			}
		}

		/** A writer thread. Takes a lock object, name descriptor, filename, and number of bytes to write */
		class Writer extends Thread
		{
			StsIOLock lock;
			private String name;
			private String filename;
			int nValues;
			FileOutputStream fos;

			Writer(StsIOLock lock, String name, String filename, int nValues)
			{
				this.lock = lock;
				setName(name);
				this.filename = filename;
				this.nValues = nValues;
				try
				{
					fos = new FileOutputStream(filename, true);
				}
				catch(Exception e) { }
			}

			public void run()
			{
				System.out.println("Starting " + name);
				byte[] bytes = new byte[nValues];
				for(int n = 0; n < nValues; n++)
					bytes[n] = (byte)n;
				try
				{
					lock.requestWrite("test write");
					fos.write(bytes);
					lock.writeAccomplished("test write");
					System.out.println("Stopping " + name + " wrote " + nValues + " bytes.");
					fos.close();
					Thread.currentThread().sleep(500);
				}
				catch(Exception e) {}
			}
		}
	}
}