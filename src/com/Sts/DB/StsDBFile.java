//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DB;

import com.Sts.Collaboration.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.channels.*;
import java.util.*;

public abstract class StsDBFile
 {
     private StsAbstractFile file;
     /** file with name db.projectName */
     private String filename = null;
     /** directory containing db file: projectDirectory/DBFiles/ */
     private String urlDirectory = null;
     /** project directory determined when file has been defined or selected.  Picked up by project on db load. */
     private String projectDirectory;

     //    private boolean checkSync = true;  // checks this dbFile against main DB to see if they are in sync
     private OutputStream os = null; // output file
     private BufferedOutputStream bos = null; // connected to output file

     protected ByteArrayOutputStream baos = null; // scratch output byte array containing transaction (copied to bos)
     protected DataOutputStream dos = null; // connected to baos

     private InputStream is = null; // input file
     private BufferedInputStream bis = null; // connected to input file
     private DataInputStream dis = null; // connected to buffered input file

     private long position = 0; // position in file
     private long cmdStartPosition;
     private long transactionStartPosition; // position before start of next command; used to set position when endDBCmd is found
     private long positionEndDB; // position after db transactions have been read (start of properties, view, info)
     private int nTransaction = 0;
     public byte status = CLOSED;

     private StsDBObjectTypeList currentDBClasses;

     private StsIOLock lock = new StsIOLock();
     private long lastTransactionTime;

     /** defines how transaction is saved */
     public byte transactionType = TRANSACTION_SAVE;

     private StsProgressBar progressBar = null;

     static public final int dbVersion = 1;
     static public final int dbSubVersion = 2;

     static public final byte CLOSED = 0;
     static public final byte WRITING = 1;
     static public final byte READING = 2;
     static public final String[] statusStrings = new String[] { "CLOSED", "WRITING", "READING"};

     static public final byte TRANSACTION_SAVE = 0;
     static public final byte TRANSACTION_PEER_SEND = 1;
     static public final byte TRANSACTION_BLOCKED = 2;
     static public final String[] transactionTypeStrings = new String[]{"Save", "Peer(s) send", "Blocked"};

     static private final int bufSize = 4096;

     static private java.text.DecimalFormat versionFormat = new java.text.DecimalFormat("##.##");

     static public final boolean collabDebug = Main.isCollabDebug;
    static public final boolean testDebug = false;
    static public final boolean debug = testDebug || Main.isDbDebug || Main.isDbCmdDebug || Main.isDbIODebug;
	static public final boolean debugPosition = Main.isDbPosDebug;

     public StsDBFile()
     {
         super();
     }

     protected StsDBFile(StsAbstractFile file, StsProgressBar progressBar) throws StsException, IOException
     {
         this.progressBar = progressBar;
         this.file = file;
         //        this.checkSync = checkSync;
         filename = file.getFilename();
         currentDBClasses = new StsDBObjectTypeList("currentDBClasses");
         if (Main.isJarDB)
         {
             urlDirectory = file.getURL().getPath();
             urlDirectory = urlDirectory.substring(0,urlDirectory.lastIndexOf("/"));
             //urlDirectory = System.getProperty("user.dirNo");
         }
         else
             urlDirectory = file.getURLDirectory();
     }

     public String getProjectDirectory()
     {
         File dbFilesDirectoryFile = new File(urlDirectory);
         return dbFilesDirectoryFile.getParent();
     }

     protected StsDBFile(StsProgressBar progressBar) throws StsException, IOException
     {
         this.progressBar = progressBar;
         currentDBClasses = new StsDBObjectTypeList("currentDBClasses");
     }

     protected StsDBFile(DataInputStream dis, StsProgressBar progressBar) throws StsException, IOException
     {
         this.dis = dis;
         this.progressBar = progressBar;
         currentDBClasses = new StsDBObjectTypeList("currentDBClasses");
     }

     protected abstract String getDBTypeName();

     static public StsObject getModelObject(Class c, int index)
     {
         StsModel model = StsObject.getCurrentModel();
         if (index < 0)
             return null;
         else
         {
             StsObject object;
             object = model.getStsClassObjectWithIndex(c, index);
             if (object != null)
                 return object;
             return model.getEmptyStsClassObject(c, index);
         }
     }

     public String getFilename()
     {
         return filename;
     }

     public String getURLDirectory()
     {
         return urlDirectory;
     }

     public String getURLPathname()
     {
         return file.getURLPathname();
     }

     public void blockTransactions()
     {
         transactionType = TRANSACTION_BLOCKED;
     }

     public void saveTransactions()
     {
         transactionType = TRANSACTION_SAVE;
     }

     public void peerSendTransactions()
     {
         transactionType = TRANSACTION_PEER_SEND;
     }

     public boolean transactionsBlocked() { return transactionType == TRANSACTION_BLOCKED; }

     static public String transactionTypeString(byte type)
     {
         return transactionTypeStrings[type];
     }

     public boolean closeReadOpenWriteAppend()
     {
         if (status == WRITING)
         {
             StsException.systemError("StsDBFile.CloseReadOpenWriteAppend() failed." +
                 " Status is not currently READING or CLOSED when attempting to close.");
             return false;
         }
         close();
         openWrite(true);
         return true;
     }

     public boolean closeReadOpenWritePosition()
     {
         if (status == WRITING)
         {
             StsException.systemError("StsDBFile.CloseReadOpenWriteAppend() failed." +
                 " Status is not currently READING or CLOSED when attempting to close.");
             return false;
         }
         close();
         if(positionEndDB == 0)
         {
            StsException.systemError(this, "closeReadOpenWritePosition", "Previous session was abnormally terminated.  Unable to recover window and property data.");
			openWrite(position);
         }
		 else
         openWrite(positionEndDB);
         if(Main.isDbPosDebug) debugCheckWritePosition("model after closeReadOpenWritePosition");
         return true;
     }
     public boolean openAndPositionToLastDBWrite()
     {
     	try
		 {
			 if(positionEndDB == 0)
			 {
				StsException.systemError(this, "repositionToLastDBWrite", "Archive didn't write positionEndDB.");
				return false;
			 }
			 if (status != CLOSED)
			 {
				 StsException.systemError(this, "DB should be closed at this point (typically after archiving).");
				 return false;
			 }
			 openWrite(true);
			 os = file.getOutputStream(true);
             if (os == null) return false;
			 FileChannel channel;
             channel = ((FileOutputStream) os).getChannel();
			 channel.position(positionEndDB);
             channel.truncate(positionEndDB);
			 if(Main.isDbPosDebug)
				StsException.systemDebug(this, "repositionToLastDBWrite", "File truncated to " + file.length());
         	 close();
         	 openWrite(true);
			 this.position = positionEndDB;
			 channel = ((FileOutputStream) os).getChannel();
			 channel.position(positionEndDB);
         	 if(Main.isDbPosDebug) debugCheckWritePosition("model after reposition to last dbWrite after archiving.");
         	return true;
     	}
		 catch(Exception e)
		 {
			 StsException.outputWarningException(this, "repositionToLastDBWrite", e);
			 return false;
		 }
     }

     public boolean openWriteNew(OutputStream outputStream)
     {
         if (status == WRITING)
         {
             StsException.systemError("StsDBFile.CloseReadOpenWriteAppend() failed." +
                 " Status is not currently READING or CLOSED when attempting to close.");
             return false;
         }
         close();
         openWrite(outputStream);
         return true;
     }

     public boolean close()
     {
         try
         {
             if (status == CLOSED)
                 return true;

             if (status == WRITING)
             {
                 bos.flush();
                 os.flush();
                 dos.flush();
                 bos = null;
                 os.close();
                 os = null;
                 baos.close();
                 baos = null;
                 dos.close();
                 dos = null;
                 if (debug && file != null)
                     System.out.println("FILE DEBUG: Closing writable file: " + filename + " length: " + file.length() +
                         " thread: " + Thread.currentThread().getName() +
                         " time: " + System.currentTimeMillis());
             }
             else if (status == READING && file != null)
             {
                 if (debug)
                     System.out.println("FILE DEBUG: Closing readable file: " + filename + " length: " + file.length() +
                         " thread: " + Thread.currentThread().getName() +
                         " time: " + System.currentTimeMillis());
                 bis.close();
                 bis = null;
                 is.close();
                 is = null;
                 dis.close();
                 dis = null;
             }

             status = CLOSED;
             // StsException.systemDebug(this, "close", "file: " + filename + " status: CLOSED");
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.close() failed."
                 + "Unable to close file " + filename, e, StsException.WARNING);
             return false;
         }
     }

     public StsDBObjectTypeList getCurrentDBClasses()
     {
         return currentDBClasses;
     }

     public long getLastTransactionTime()
     {
         return lastTransactionTime;
     }

     public StsDBType getCurrentDBType(Class c)
     {
         return currentDBClasses.getDBType(c);
     }

     public StsDBType getCurrentDBType(String s)
     {
         return currentDBClasses.getDBType(s, null);
     }

     public StsDBType getCurrentDBType(Object object)
     {
         return currentDBClasses.getDBType(object);
     }
/*
     public void setDBTypesTemporary(boolean temporary)
     {
         StsDBType.setTypeIsTemporary(temporary);
     }
*/
     public void removeTemporaryDBTypes()
     {
         currentDBClasses.removeTemporaryTypes();
     }

     public boolean commitCmdList(String name, ArrayList cmdList)
     {
         if (cmdList == null)
             return false;
         if (status != WRITING)
             return false;
         if (transactionType == TRANSACTION_BLOCKED)
             return false;

         try
         {
             ByteArrayOutputStream outputBytes = new ByteArrayOutputStream();
             StsDBOutputStream dbOutputStream = new StsDBOutputStream(this, outputBytes, currentDBClasses);

             // write out commands to dbOutputStream
             int nCmds = cmdList.size();
             for (int i = 0; i < nCmds; i++)
             {
                 StsDBCommand cmd = (StsDBCommand) cmdList.get(i);
                 write(dbOutputStream, cmd);
             }

             writeTransaction(outputBytes, dbOutputStream, name);

             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.commitCmdList() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     public boolean commitCmd(String name, StsDBCommand cmd)
     {
         if (status != WRITING)
             return false;
         if (transactionType == TRANSACTION_BLOCKED)
             return false;
         if (cmd == null)
             return false;

         try
         {
             ByteArrayOutputStream outputBytes = new ByteArrayOutputStream();
             StsDBOutputStream dbOutputStream = new StsDBOutputStream(this, outputBytes, currentDBClasses);
			 transactionStartPosition = position;
			 if(debugPosition)
			 {
				 StsException.systemDebug(this, "commitCmd", "Transaction: " + nTransaction + " transaction start set: " + position);

			 }
             write(dbOutputStream, cmd);

             writeTransaction(outputBytes, dbOutputStream, name);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.commitCmdList() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     /** Used by collaboration to read bytes received in transactionMessage from leader */
     public boolean readTransaction(byte[] transactionBytes, String description, byte incomingTransactionType, StsModel model)
     {
         try
         {
             ByteArrayInputStream bais = new ByteArrayInputStream(transactionBytes);
             DataInputStream dis = new DataInputStream(bais);
             transactionType = TRANSACTION_BLOCKED; // block any  writes to db in read process (which can possibly create db commands)
             boolean readOk = readObjects(dis);
             bais.close();
             dis.close();
             return readOk;
         }
         catch (EOFException eof)
         {
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.readObjects() failed for transaction " + description, e, StsException.WARNING);
             return false;
         }
         finally
         {
             model.refreshObjectPanel();
             try
             {
                 if (transactionType == TRANSACTION_SAVE)
                 {
                     bos.write(transactionBytes);
                     if (debug)
                         System.out.println("FILE DEBUG: for transaction " + description + " readTransaction wrote " + transactionBytes.length + " bytes to file " +
                             filename);
                 }
             }
             catch (Exception e)
             {
                 StsException.systemError("StsDBFile.readTransaction() failed for transaction " + description + " writing to db file.");
             }
         }
     }

     synchronized public boolean copyTo(StsDBFile newDb)
     {
         try
         {
             if (newDb == null)
                 return false;
             if (!openReadAndCheckFile(false))
                 return false;
             if (!newDb.openWrite())
                 return false;

             byte[] buffer = new byte[bufSize];

             int bytesRead;
             while ((bytesRead = dis.read(buffer)) != -1)
                 newDb.write(buffer, 0, bytesRead);

             close();
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.copyTo() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     /**
      * used for collaboration when new peer reads the inputStream from the leader
      * and writes the data to a file which is subsequently opened and read to start.
      */
     static public StsFile copyToFile(DataInputStream dataInputStream)
     {
         try
         {
             String pathname = "C:/temp/db.temp";
             File f = new File(pathname);
             if (f.exists())
                 f.delete();
             //			File f = File.createTempFile("db.", ".temp");
             //			String pathname = f.getPath();
             StsFile file = StsFile.constructor(pathname);
             if (collabDebug)
                 System.out.println("Peer is using temp file: " + pathname);

             OutputStream os = file.getOutputStream(false);
             if (os == null)
                 return null;
             BufferedOutputStream bos = new BufferedOutputStream(os, bufSize);

             long nBytesSent = dataInputStream.readLong();
             if (collabDebug)
                 System.out.println("nBytes sent by leader: " + nBytesSent);

             byte[] buffer = new byte[bufSize];
             long nTotalBytesRead = 0;
             int bytesRead;
             while ((bytesRead = dataInputStream.read(buffer)) != -1)
             {
                 nTotalBytesRead += bytesRead;
                 // System.out.println("bytes received: " + bytesRead + " total bytes so far: " + nTotalBytesRead);
                 bos.write(buffer, 0, bytesRead);
                 if (nTotalBytesRead >= nBytesSent)
                     break;
             }
             if (collabDebug)
                 System.out.println("total bytes received: " + nTotalBytesRead);
             if (debug)
                 System.out.println("FILE DEBUG: copyToFile copied " + nTotalBytesRead + " bytes to file " + pathname);

             bos.flush();
             os.close();
             return file;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.copyTo() failed.", e, StsException.WARNING);
             return null;
         }
     }

     /**
      * Used in collaboration when the leader wants to write the database to an outputStream
      * which is on one end of a socket.  Peer is on the other and uses copyToFile() to read inputStream
      */
     public boolean writeToOutputStream(OutputStream os)
     {
         try
         {
             close();
             long fileSize = file.length();
             if (collabDebug)
                 System.out.println("File size transferred by leader: " + fileSize);
             byte[] buffer = new byte[bufSize];
             DataOutputStream dos = new DataOutputStream(os);
             long nTotalBytesWritten = 0;
             try
             {
                 dos.writeLong(fileSize);
                 lock.requestRead("writeToOutputStream(os)");
                 openRead();
                 int lastBytesRead = 0, bytesRead;
                 while ((bytesRead = dis.read(buffer)) != -1)
                 {
                     dos.write(buffer, 0, bytesRead);
                     nTotalBytesWritten += bytesRead;
                     if (collabDebug)
                         System.out.println("bytes written: " + bytesRead + " total bytes written so far: " + nTotalBytesWritten);
                     lastBytesRead = bytesRead;
                 }
                 if (collabDebug)
                     System.out.println("Total bytes written by leader: " + nTotalBytesWritten);
                 return true;
             }
             catch (Exception e)
             {
                 StsException.outputException("StsDBFile.writeToOutputStream() failed.",
                     e, StsException.WARNING);
                 return false;
             }
             finally
             {
                 lock.readAccomplished("writeToOutputStream(os)");
                 closeReadOpenWriteAppend();
             }
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.writeToOutputStream() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     /**
      * Used in collaboration when the leader wants to write the database to an outputStream
      * which is on one end of a socket.  Peer is on the other and uses copyToFile() to read inputStream
      */
     public boolean writeBytesToOutputStream(byte[] outputBytes, OutputStream os)
     {
         try
         {
             DataOutputStream dos = new DataOutputStream(os);
             int size = outputBytes.length;
             dos.writeInt(size);
             dos.writeLong(size);
             dos.write(outputBytes, 0, size);
             dos.flush();
             System.out.println("Total transaction bytes written by leader: " + size);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.writeToOutputStream() failed.",
                 e, StsException.WARNING);
             return false;
         }
     }

     /** open file for reading from the beginning of the file */
     private boolean openRead()
     {
         try
         {
             if (debug)
                 System.out.println("FILE DEBUG: Opening file to read: " + filename + " length: " + file.length() +
                     " thread: " + Thread.currentThread().getName() +
                     " time: " + System.currentTimeMillis());
             is = file.getInputStream();
             bis = new BufferedInputStream(is, bufSize);
             dis = new DataInputStream(bis);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openReadAndCheck() failed.", e, StsException.WARNING);
             return false;
         }
     }

     private boolean openRead(InputStream is)
     {
         try
         {
             if (debug)
                 System.out.println("FILE DEBUG: Opening input stream to read " +
                     " thread: " + Thread.currentThread().getName() +
                     " time: " + System.currentTimeMillis());
             bis = new BufferedInputStream(is, bufSize);
             dis = new DataInputStream(bis);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openReadAndCheck() failed.", e, StsException.WARNING);
             return false;
         }
     }
     /** open file for reading from the beginning of the file */
 /*
     private boolean openRead(InputStream inputStream)
     {
         try
         {
             is = inputStream;
             ((FileInputStream) is).getChannel().position(0);
             bis = new BufferedInputStream(is, bufSize);
             dis = new DataInputStream(bis);
             if (!isDBFile(dis)) return false;
             status = READING;
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openReadAndCheck() failed.", e, StsException.WARNING);
             return false;
         }
     }
*/
     /* open file for reading, check that it begins with correct string ID, leaving file position after ID. */
     protected boolean openReadAndCheckFile(boolean check)
     {
         try
         {
             if (!openRead())
                 return false;
             if (!isDBFile(dis, check))
                 return false;
             status = READING;
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openReadAndCheckFile() failed.", e, StsException.WARNING);
             return false;
         }
     }

     protected boolean openWrite()
     {
         return openWrite(false);
     }

     private boolean openWrite(boolean append)
     {
         try
         {
             if (debug) System.out.println("FILE DEBUG: Opening file to write: " + filename +
                 " append: " + append + " length: " + file.length() +
                 " thread: " + Thread.currentThread().getName() +
                 " time: " + System.currentTimeMillis());
             os = file.getOutputStream(append);
             if (os == null) return false;
             bos = new BufferedOutputStream(os, bufSize);

             baos = new ByteArrayOutputStream(bufSize); // scratch array holding transaction
             dos = new DataOutputStream(baos);

             if (!append)
             {
                 dos.writeUTF(getDBTypeName() + " " + StsModel.version);
                 dos.writeInt(dbVersion);
                 dos.writeInt(dbSubVersion);
                 int datalen = baos.size() + 4;
                 dos.writeInt(datalen);
                 dos.flush();
                 outputBytes("file header");
                 position = datalen;
                 if (debugPosition)
                 {
                     debugCheckWritePosition("openWrite(append = false). DB File header size: " + datalen + " for file " + filename);
                 }
             }
             status = WRITING;

             return true;
         }
         catch (FileNotFoundException fnfe)
         {
             String errorMessage = "Failed to find or access database file " + filename;
             new StsMessage(null, StsMessage.ERROR, errorMessage);
             return false;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openWrite() failed.", e, StsException.WARNING);
             return false;
         }
     }

     private boolean openWrite(long position)
     {
         if (position == 0) return openWrite();
         try
         {
             this.position = position;
             os = file.getOutputStream(true);
             if (os == null) return false;

             FileChannel channel = ((FileOutputStream) os).getChannel();
             channel.truncate(position);

             if (debugPosition)
             {
                 debugCheckWritePosition("openWrite(position)");
                 System.out.println("FILE DEBUG: Open output file and position " + filename +
                     " position: " + position + " length: " + file.length() +
                     " thread: " + Thread.currentThread().getName() +
                     " time: " + System.currentTimeMillis());
             }
             
             bos = new BufferedOutputStream(os, bufSize);

             baos = new ByteArrayOutputStream(bufSize); // scratch array holding transaction
             dos = new DataOutputStream(baos);

             status = WRITING;

             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openWrite() failed.", e, StsException.WARNING);
             return false;
         }
     }

     private boolean openWrite(OutputStream outputStream)
     {
         try
         {
             if (outputStream == null) return false;
             os = outputStream;
             bos = new BufferedOutputStream(os, bufSize);

             baos = new ByteArrayOutputStream(bufSize); // scratch array holding transaction
             dos = new DataOutputStream(baos);

             dos.writeUTF(getDBTypeName() + " " + StsModel.version);
             dos.writeInt(dbVersion);
             dos.writeInt(dbSubVersion);
             int datalen = baos.size() + 4;
             dos.writeInt(datalen);
             dos.flush();
             outputBytes("file header");
             if (debug) System.out.println("DB File header size: " + datalen + " for file " + filename);
             status = WRITING;
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.openWrite() failed.", e, StsException.WARNING);
             return false;
         }
     }

     private void writeTransaction(ByteArrayOutputStream outputBytes, StsDBOutputStream dbOutputStream, String name)
     {
         if (status != WRITING)
             return;
         if (transactionType == TRANSACTION_BLOCKED)
             return;

         try
         {
             lock.requestWrite("writeTransaction(" + name + ")");
             if (debug) debugStartWriteTransaction(name);
             // Write out the transaction begin header first
             ByteArrayOutputStream headerBytes = new ByteArrayOutputStream();
             StsDBOutputStream headerOutputStream = new StsDBOutputStream(this, headerBytes, currentDBClasses);

             write(headerOutputStream, new StsStartTransactionCmd(name));
             headerOutputStream.flush();

             int datalen = headerBytes.size();
             if (debug) System.out.println("Writing out " + datalen + " header bytes to " + filename);
             dos.writeInt(datalen);
             byte[] bytesToWrite = headerBytes.toByteArray();
             dos.write(bytesToWrite);
             position += (4 + datalen);
             dos.flush();
             outputBytes(new String("Transaction header " + nTransaction + ": " + name));
             if (debugPosition) debugCheckWritePosition("wrote transaction header");

             // Writing commands may have caused some requests to save dbClasses. Write these second
             ByteArrayOutputStream classDefinitionBytes = new ByteArrayOutputStream();
             StsDBOutputStream dbClassDefinitionOutputStream = new StsDBOutputStream(this, classDefinitionBytes, currentDBClasses);
             if (debug) currentDBClasses.debugPrintNewClasses();
             writeDBClasses(dbClassDefinitionOutputStream, currentDBClasses.getNewClasses());
             currentDBClasses.clearNewClasses();
             dbClassDefinitionOutputStream.flush();
             datalen = classDefinitionBytes.size();
             if (debug) System.out.println("Writing out " + datalen + " class definition bytes to " + filename);
             dos.writeInt(datalen);
             bytesToWrite = classDefinitionBytes.toByteArray();
             dos.write(bytesToWrite);
             dos.flush();
             outputBytes(new String("Transaction Class Definitions " + nTransaction + ": " + name));
             position += (4 + datalen);
             if (Main.isDbIODebug) debugCheckWritePosition("writeTransaction class definitions");

             // now commit the commands to disk
             dbOutputStream.flush();
             datalen = outputBytes.size();
             if (debug) System.out.println("Writing out " + datalen + " data bytes to " + filename);
             dos.writeInt(datalen);
             bytesToWrite = outputBytes.toByteArray();
             dos.write(bytesToWrite);
             dos.flush();
             outputBytes(new String("Transaction Objects" + nTransaction + ": " + name));
             position += (4 + datalen);
             if (Main.isDbIODebug) debugCheckWritePosition("writeTransaction objects");

             if (debug) debugEndWriteTransaction(name);
             nTransaction++;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.writeTransaction() failed.",
                 e, StsException.WARNING);
         }
         finally
         {
             lock.writeAccomplished("writeTransaction(" + name + ")");
			 removeTemporaryDBTypes();
         }
     }

     public void debugCheckWritePosition(String message)
     {
         System.out.println(message + " position: " + position);
         // long filePosition = getWriteFilePosition();
         //System.out.println(message + " filePosition: " + filePosition + " position: " + position);
         /*
            if(filePosition != position)
             {
                 StsException.systemError(this, "debugCheckWritePosition", "filePosition and position don't agree.");
             }
         */
     }

     public void debugCheckReadPosition(String message)
     {
         System.out.println(message + " position: " + position);
         // long filePosition = getReadFilePosition();
         // System.out.println(message + " filePosition: " + filePosition + " position: " + position);
         /*
        if(filePosition != position)
         {
             StsException.systemError(this, "debugCheckReadPosition", "filePosition and position don't agree.");
         }
         */
     }

     private long getWriteFilePosition()
     {
         try
         {
             return ((FileOutputStream) os).getChannel().position();
         }
         catch (Exception e)
         {
             StsException.systemError(this, "getWriteFilePosition");
             return 0;
         }
     }

     private long getReadFilePosition()
     {
         try
         {
             return ((FileInputStream) is).getChannel().position();
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "getReadFilePosition", e);
             return 0;
         }
     }

     public void debugStartWriteTransaction(String name)
     {
         debugCheckWritePosition("START WRITE TRANSACTION " + nTransaction + ": " + name);

     }

     public void debugEndWriteTransaction(String name)
     {
         debugCheckWritePosition("END WRITE TRANSACTION " + nTransaction + ": " + name);
     }

     /** debug display filePosition;  position is not in sync because number of bytes read for this transaction are added yet */
     public void debugStartReadTransaction(String name)
     {
         debugCheckReadPosition("START READ TRANSACTION " + nTransaction + ": " + name);
     }

     public void debugEndReadTransaction(String name)
     {
         debugCheckReadPosition("END READ TRANSACTION " + nTransaction + ": " + name);
     }

     private void outputBytes(String description)
     {
         try
         {
             byte[] bytes = baos.toByteArray();
             if (bytes == null) return;
             if (bos != null && transactionType == TRANSACTION_SAVE)
             {
                 if (debug)
                     System.out.println("FILE DEBUG: For " + description + " wrote " + bytes.length + " bytes to file " + filename);
                 bos.write(bytes);
                 bos.flush();
                 os.flush();
             }
             baos.reset();

             StsCollaboration collaboration = StsCollaboration.getCollaboration();
             if (collaboration != null)
             {
                 if (collaboration.isLeader())
                 {
                     collaboration.broadcastTransaction(bytes, description, transactionType);
                     if (debug)
                         System.out.println("                 sending to peer(s)");
                 }
             }
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.outputBytes() failed writing " + description,
                 e, StsException.WARNING);
             return;
         }
     }

     synchronized private boolean write(StsDBOutputStream dbOutputStream, StsDBCommand cmd)
     {
         try
         {
             if (cmd == null)
                 return false;
             try
             {
                 //				byte index = cmd.getDBCommandClassIndex();
                 //				dbOutputStream.writeByte(index);
                 // cmdStartPosition = position;
                 cmd.write(dbOutputStream);
             }
             catch (Exception e)
             {
                 StsException.outputException("StsDBFile.write(StsDBOutputStream, StsDBCommand) failed.",
                     e, StsException.WARNING);
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.write(obj) failed.", e, StsException.WARNING);
             return false;
         }
     }

     private void readDBCommand(StsDBInputStream dbInputStream) throws IOException, InstantiationException, IllegalAccessException
     {
         // cmdStartPosition = position;
         byte dbCommandClassIndex = dbInputStream.readByte();
         Class dbCommandClass = StsDBCommand.getDBCommandClass(dbCommandClassIndex);
         if (dbCommandClass == null) return;
         StsDBCommand command = (StsDBCommand) dbCommandClass.newInstance();
         command.read(dbInputStream);
     }
/*
     private void writeDBClasses(StsDBOutputStream dbClassDefinitionOutputStream, Object[] dbClasses)
     {
         for (int n = 0; n < dbClasses.length; n++)
         {
             StsDBType dbClass = (StsDBType) dbClasses[n];
             if (dbClass.representsStsSerialisableObject())
             {
                 StsDBTypeClass dbc = (StsDBTypeClass) dbClass;
                 if (!dbc.getMatchesOnDiskVersion())
                 {
                     if (debug)
                         System.out.println(this.getFilename() + " - Writing " + dbc.getTypeName());

                     StsDBTypeCmd cmd = new StsDBTypeCmd(dbc);
                     write(dbClassDefinitionOutputStream, cmd);
                 }
             }
         }
	}
*/
     private void writeDBClasses(StsDBOutputStream dbClassDefinitionOutputStream, Object[] dbClasses)
     {
         for (int n = 0; n < dbClasses.length; n++)
         {
             StsDBType dbClass = (StsDBType) dbClasses[n];
             if (dbClass.representsStsSerialisableObject())
             {
                 StsDBTypeClass dbc = (StsDBTypeClass) dbClass;
                 boolean temporary = dbc.isTemporary();
                 boolean matchesDisk = dbc.getMatchesOnDiskVersion();
                 if (!matchesDisk || temporary)
                 {
                     if (debug) System.out.println(this.getFilename() + " - Writing " + dbc.getTypeName());
                     StsDBTypeCmd cmd = new StsDBTypeCmd(dbc);
                     write(dbClassDefinitionOutputStream, cmd);
                     dbc.setMatchesOnDiskVersion(true);
                 }
             }
         }
     }
     protected boolean readObjects()
     {
         String lockString = "readFileObjects()";
         try
         {
             if (status != READING) return false;
             lock.requestRead(lockString);
             return readObjects(dis);
         }
         catch (EOFException eof)
         {
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsDBFile.readObjects() failed.", e, StsException.WARNING);
             return false;
         }
         finally
         {
             lock.readAccomplished(lockString);
             if (progressBar != null)
             {
                 progressBar.setValueAndStringImmediate(100, "Completed loading database");
                 StsSplashScreen.instance.paintImmediately();
             }
         }
     }

     static final int nProgressIntervals = 100;

     synchronized private boolean readObjects(DataInputStream dis) throws EOFException, IOException
     {
         int datalen = 0;
         StsDBInputStream istrm = null;

         try
         {
             double totalBytes = dis.available();
             int progressPositionInc = (int)(totalBytes/nProgressIntervals);
             long nextProgressPosition = position;

             int byteReadInc = 0;
             if (progressBar != null)
             {
                 progressBar.initializeImmediate(100);
                 progressBar.setStringImmediate("Reading " + this.filename);
                 StsSplashScreen.instance.paintImmediately();
             }
             boolean readOk = true;
             while (dis.available() > 0)
             {
                 // get the length of the data chunk
                 datalen = dis.readInt();
                 if (debugPosition)
					 debugCheckReadPosition("Reading in data chunk " + (datalen + 4) + " bytes from buffer");

                 if (progressBar != null && position >= nextProgressPosition)
                 {
                     progressBar.setValueImmediate((int) (100 * position / totalBytes));
                     nextProgressPosition = position + progressPositionInc;
                     StsSplashScreen.instance.paintImmediately();
                 }

                 // retrieve the object data and store it
                 byte[] data = new byte[datalen];
                 dis.readFully(data);
                 // position += (4 + datalen);
                 if (debugPosition)
					 debugCheckReadPosition("Create object(s) byte stream of length " + (datalen + 4));
                 // create a byte stream to contain the input object
                 ByteArrayInputStream bytes = new ByteArrayInputStream(data);

                 // create a stream to translate the object
                 istrm = getDBInputStream(bytes, istrm);
                 try
                 {
                     while (istrm.available() > 0)
                     {
                         readDBCommand(istrm);
                     }
					 position += (4 + datalen);
                 }
                 catch (EOFException eofe)
                 {
                     StsException.systemError(this, "readObjects", " EOF unexpectedly encountered.");
                     return false;
                 }
                 catch (Exception e)
                 {
                     StsException.outputWarningException(this, "readObjects", " Attempting to continue with next data block.", e);
                     readOk = false;
                 }
                 //                position += (4 + datalen); // 4 is for datalen int
             }
             // Now need to set the currentDBClasses to have the sames indexes as the inputDBClasses
             currentDBClasses = istrm.getInputTypeList();
             currentDBClasses.removeTemporaryTypes();
             if(debug)
             {
                 System.out.println("DBTypes after read completed:");
                 currentDBClasses.printDBTypes();
             }
             return readOk;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "readObjects", "failed reading from input: " + "Couldn't read " + datalen + " bytes.", e);
             return false;
         }
         finally
         {
             istrm.initializeObjects(StsModel.getCurrentModel());
         }
     }

     public void setPositionEndDB()
     {
         if(Main.isDbPosDebug) StsException.systemDebug(this, "setPositionEndDB", "positionEndDB set to start of transaction: " + transactionStartPosition);
         positionEndDB = transactionStartPosition;
         if(positionEndDB == 0)
            StsException.systemError(this, "setPositionEndDB", "positionEndDB should not be zero!");
     }

     protected StsDBInputStream getDBInputStream(InputStream in, StsDBInputStream oldDBInputStream) throws IOException
     {
         return new StsDBInputStream(this, in, oldDBInputStream, currentDBClasses);
     }

     private boolean isDBFile(DataInputStream dbInputStream, boolean check)
     {
         String inDBTypeName = null;
         int inDBVersion;
         int inDBSubVersion;
         try
         {
             inDBTypeName = dbInputStream.readUTF();
             inDBVersion = dbInputStream.readInt();  // not currently used
             inDBSubVersion = dbInputStream.readInt(); // not currently used
             int headerSize = dbInputStream.readInt();
             if(check && !dbTypeNameOk(inDBTypeName)) return false;

             position = headerSize;
			 if(debugPosition)
				 debugCheckReadPosition("Checking file. headerSize: " + headerSize);
             if (debug)
             {
                 System.out.print("Reading database of type " + inDBTypeName + " Version " + formatVersion(inDBVersion, inDBSubVersion));
                 System.out.println(", Current Version is " + formatVersion(dbVersion, dbSubVersion));
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.systemError("Couldn't read file type. This is not a db file.");
             return false;
         }
     }

     private boolean dbTypeNameOk(String inDbTypeString)
     {
        String[] tokens = StsToolkit.getTokens(inDbTypeString);
        int nTokens = tokens.length;
        String inDBTypeName = tokens[0];
        if (!inDBTypeName.equalsIgnoreCase(getDBTypeName()))
            return false;
        if(nTokens < 2)  return true;
        String dbModelVersion = tokens[1];
        System.out.println("Reading database built with version " + dbModelVersion);

        if(StsStringUtils.compareTokenStrings(StsModel.version, dbModelVersion, ".") >= 0) return true;

        return StsYesNoDialog.questionValue(null, "This S2S model version: " + StsModel.version + "\n" +
            " is older than the version which built this db: " + dbModelVersion + ".\n" +
            " Do you wish to continue?");
     }


     private String formatVersion(float version, float subVersion)
     {
         return versionFormat.format(version + (subVersion / 100f));
     }

     private void write(byte[] bytes, int offset, int length) throws IOException
     {
         bos.write(bytes, offset, length);
     }

     public void setLastTransactionTime(long lastTransactionTime)
     {
         this.lastTransactionTime = lastTransactionTime;
     }

     public long getPosition()
     {
         return position;
     }

     public void setStartTransactionPosition()
     {
         if(nTransaction != 0 && position == 0)
             StsException.systemDebug(this, "setStartTransactionPosition", "For nTransaction " + nTransaction + " transactionStartPosition should not be zero!");
         transactionStartPosition = position;
		 if(debugPosition) StsException.systemDebug(this, "setStartTransactionPosition", "Transaction: " + nTransaction + " position " + position);
     }

     public void incrementNTransaction()
     {
         nTransaction++;
     }
 }