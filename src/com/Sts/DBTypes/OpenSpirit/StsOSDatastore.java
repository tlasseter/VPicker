package com.Sts.DBTypes.OpenSpirit;

import com.Sts.Actions.Import.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.DataTransfer.*;
import com.Sts.UI.Progress.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 10, 2008
 * Time: 7:37:23 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsOSDatastore extends StsObject implements StsDatastoreFace
{
    boolean debug = true;
    String datastoreName = null;
    String projectName = null;
    boolean projectLimit = true;
    StsOpenSpiritImport ospImport = null;

    abstract public Object[] getProjectObjects(StsProgressPanel panel, String[] wellboreList);

    /**
     * Default empty constructor
     */
    public StsOSDatastore()
    {
    }

    /**
     * Constructor for finding available datastores from OpenSpirit. The
     * datastore will be added later.
     * @param ospImport is the input OpenSpirit instance.
     */
    public StsOSDatastore(StsOpenSpiritImport ospImport)
    {
    	this.ospImport = ospImport;

    	// connect to the OpenSpirit instance if necessary
    	if (this.ospImport == null)
    	{
    		this.ospImport = new StsOpenSpiritImport(currentModel);
    		this.ospImport.initializeOpenSpirit();
    	}
    }

    /**
     * Constructor when the datastoreName to use is already known.
     * @param ospImport is the input OpenSpirit instance.
     * @param datastoreName is the OpenSpirit datastore name for this
     * StsOSDatastore instance.
     */
    public StsOSDatastore(StsOpenSpiritImport ospImport, String datastoreName)
    {
    	this.ospImport = ospImport;
        this.datastoreName = datastoreName;
        System.out.println("StsOSDatastore construction: Set datastore name to " + this.datastoreName);
    	// connect to the OpenSpirit instance if necessary
    	if (this.ospImport == null)
    		this.ospImport = new StsOpenSpiritImport(currentModel);
    	if (!ospImport.getOpenSpiritInstance().isConnected())
    	{
    		ospImport.initializeOpenSpirit();
    	}

    }

    public boolean initialize(StsModel model)
    {
        return initialize();
    }

    public boolean initialize()
    {
        return true;
    }


    public String[] getDatastoreNames() { return new String[] { datastoreName }; }
    public void setDatastoreNames(String[] names) { datastoreName = names[0]; }

    public void setProjectLimit(boolean limit)
    {
        projectLimit = limit;
    }

    public String getProjectName()
    {
    	return projectName;
    }

    public void setOpenSpiritImport(StsOpenSpiritImport osp)
    {
    	ospImport = osp;
    }

    public StsOpenSpiritImport getOpenSpiritImport()
    {
    	return ospImport;
    }

 	public String getDatastoreName()
	{
		return datastoreName;
	}

	public void setDatastoreName(String datastoreName)
	{
		this.datastoreName = datastoreName;
        dbFieldChanged("datastoreName", datastoreName);
    }

	/**
 	 * Get the project names from the datastore as an array of
 	 * datastoreName - projectName Strings.
 	 */
    public String[] getProjectNames()
    {
    	if (ospImport.getOpenSpiritInstance().isConnected() &&
    			datastoreName != null && datastoreName.length() > 0)
    	{
    		String[] ds = new String[1];
    		ds[0] = datastoreName;
    		return ospImport.getProjectNames(ds);
    	}
    	else
    		return new String[0];
    }

    /**
     * Set the project name and QueryScope from the input
     * datastore - projectName string
     */
    public void setProjectName(String datastoreProjectName)
    {
    	String[] name = new String[1];
    	int iDash = datastoreProjectName.indexOf("-");
    	name[0] = datastoreProjectName.substring(iDash+2);
    	projectName = name[0];
        datastoreName = datastoreProjectName.substring(0,iDash-1);
    	name[0] = datastoreProjectName;
        //if(debug)
            System.out.println("setProjectName:Setting query scope to: " + name[0]);
    	ospImport.setQueryScope(name);
        if(debug) System.out.println("Successfully setQueryScope for project.");
        dbFieldChanged("projectName", projectName);
    }

    public String toString() { return datastoreName; }
}
