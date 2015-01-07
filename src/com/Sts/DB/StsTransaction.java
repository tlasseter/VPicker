
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DB;

import com.Sts.DB.DBCommand.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsTransaction
{
//    static private StsTransaction currentTransaction = null;
//    static private int nTransaction = 0;

    private String name = null;
 	private ArrayList cmdList = null;

    public StsTransaction(String name)
    {
//        if(Main.isDbDebug) System.out.println("TRANSACTION: " + name);
    	this.name = name;
    	cmdList = new ArrayList(10);
//        currentTransaction = this;
//        nTransaction++;
    }

    public void add(StsDBCommand cmd)
    {
        if(cmdList == null)
        {
            StsException.systemError("StsTransaction.add(cmd) failed, cmdList is null.");
            return;
        }
		cmdList.add(cmd);
    }

    public String getName()
    {
    	return name;
    }

    // not currently used, but would be nice though aborting a single command
    // can make model inconsistent
   	public boolean abortLastCmd() throws StsException
    {
        int nCmds = cmdList.size();
        if(nCmds == 0) return false;
        StsDBCommand cmd = (StsDBCommand) cmdList.get(nCmds-1);
        cmd.abort();
        return true;
    }

   	public boolean abort() throws StsException
    {
        if(cmdList == null) return false;
        int nCmds = cmdList.size();
        if(nCmds == 0) return false;
		for(int n = nCmds; n > 0; n--)
		{
			StsDBCommand cmd = (StsDBCommand) cmdList.get(n-1);
			cmd.abort();
		}
        cmdList.clear();
//        cmdList = null;
        return true;
    }

    public boolean commit(StsModel model)
    {
        if(cmdList == null || cmdList.size() == 0) return false;

        try
        {
            StsDBFile db = model.getDatabase();
            if(db == null) return false;
            db.commitCmdList(name, cmdList);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsTransaction.commit() failed.",
                e, StsException.WARNING);
            return false;
        }
        finally
        {
            cmdList = null;
//            currentTransaction = null;
        }
    }

    public void clear()
    {
        if(cmdList != null) cmdList.clear();
    }

/*
    static public void currentCommit(String name, StsModel model) throws StsException
    {
		getCurrent().commit(model);
	}
*/
/*
    static public void currentCommit(String name) throws StsException
    {
		getCurrent().commit();
	}
*/
/*
    static public void currentAbort()
    {
        try
        {
		    getCurrent().abort();
        }
        catch(Exception e)
        {
            StsException.outputException("StsTransaction.currentAbort() failed.",
                                         e, StsException.WARNING);
        }
    }
*/
    public boolean isEmpty()
    {
    	return cmdList.isEmpty();
    }
}


