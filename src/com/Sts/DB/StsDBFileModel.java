package com.Sts.DB;

import com.Sts.DB.DBCommand.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsDBFileModel extends StsDBFile
{
    static public final boolean debug = false;
    public StsDBFileModel()
	{
		super();
	}

	private StsDBFileModel(StsAbstractFile file, StsProgressBar progressBar) throws StsException, IOException
	{
		super(file, progressBar);
	}

	private StsDBFileModel(DataInputStream dis, StsProgressBar progressBar) throws StsException, IOException
	{
		super(dis, progressBar);
	}

	static public boolean fileOK(StsAbstractFile file, StsProgressBar progressBar)
	{
		StsDBFileModel dbFile = openRead(file, progressBar, true);
		if (dbFile == null) return false;
		dbFile.close();
		return true;
	}

	static public StsDBFileModel openWrite(StsAbstractFile file, StsProgressBar progressBar)
	{
		try
		{
			StsDBFileModel dbFile = new StsDBFileModel(file, progressBar);
			if (!dbFile.openWrite())
			   return null;
			return dbFile;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFile.openWrite() failed. File: " + file.getFilename());
			return null;
		}
	}

	static public StsDBFileModel openRead(StsAbstractFile file, StsProgressBar progressBar, boolean check)
	{
		try
		{
			StsDBFileModel dbFile = new StsDBFileModel(file, progressBar);
			if (! dbFile.openReadAndCheckFile(check))
			{
				return null;
			}
			return dbFile;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFile.openReadAndCheck() failed. File: " + file.getFilename());
			return null;
		}
	}

	static public StsDBFileModel openRead(DataInputStream dis, StsProgressBar progressBar)
	{
		try
		{
			return new StsDBFileModel(dis, progressBar);
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFile.openReadAndCheck(dis) failed.");
			return null;
		}
	}

	protected String getDBTypeName()
	{
		return "S2S-DB-MODEL";
	}

    public synchronized boolean writeModel(StsModel model)
	{
		StsSaveModelCmd cmd = new StsSaveModelCmd(model);
		return commitCmd("Save Model Transaction", cmd);
	}

	public synchronized boolean readModel()
	{
        if(debug) System.out.println("Reading model DB");
        boolean readOk = readObjects();
//        removeTemporaryDBTypes();
        return readOk;
	}

    static public void main(String[] args)
    {
        testArrayChangeCmd("c:\\stsdev\\c76\\arrayChgCmd");
    }
    static private void testArrayChangeCmd(String pathname)
    {
        try
        {
            Main.isDbDebug = true;
            Main.usageTracking = false;
            StsModel model = new StsModel();
            StsObject.setCurrentModel(model);
            StsFile file = StsFile.constructor(pathname);
            StsDBFileModel dbFile = new StsDBFileModel(file, null);
            model.setDatabase(dbFile);
            StsVelocityProfile profile = new StsVelocityProfile(false);
            profile.setProfilePoints(new StsPoint[] { new StsPoint(0, 0), new StsPoint(1, 1) });
            profile.addToModel();
            StsPoint newPoint = new StsPoint(2, 2);
            StsArrayInsertCmd cmd = new StsArrayInsertCmd(profile, newPoint, "profilePoints", 0, false);
            model.addTransactionCmd("arrayChange", cmd);

            model.close();

            model = new StsModel();
            StsObject.setCurrentModel(model);
            dbFile = new StsDBFileModel(file, null);
            dbFile.readModel();

            System.out.println("");
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}
