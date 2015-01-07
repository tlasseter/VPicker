package com.Sts.DBTypes.OpenSpirit;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.OSWell.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 10, 2008
 * Time: 7:49:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsOSWellDatastore extends StsOSDatastore
{
    static StsOSWellDatastore datastore = null;

/*    private StsOSWellDatastore()
    {
        super("OpenSpiritDatastore");
        int nWells = 5;
        wells = new StsOSWell[nWells];
        for(int n = 0; n < nWells; n++)
            wells[n] = new StsOSWell("OSwell-" + n);
    }
*/
    /**
     * Default empty constructor
     */
    public StsOSWellDatastore()
    {
    }

    public StsOSWellDatastore(StsOpenSpiritImport ospImport)
    {
    	super(ospImport);
    }

    public StsOSWellDatastore(StsOpenSpiritImport ospImport, String dataStoreName)
    {
		super(ospImport, dataStoreName);
    }

    static public StsOSWellDatastore getInstance()
    {
        return datastore;
    }

    public Object[] getObjects()
    {
        return getProjectObjects(null, null);
    }

    /**
     * Set the project name and QueryScope from the input
     * datastore - projectName string
     */
    public void setProjectName(String datastoreProjectName)
    {
        super.setProjectName(datastoreProjectName);

        // Run get objects

    }

    public void runGetProjectObjects(StsOSWellObjectTransferPanel _transferPanel, StsProgressPanel _progressPanel, String[] _wellboreList)
    {
        final StsOSWellObjectTransferPanel transferPanel = _transferPanel;
        final StsProgressPanel progressPanel = _progressPanel;
        final String[] wellboreList = _wellboreList;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                Object[] objects = getProjectObjects(progressPanel, wellboreList);
                transferPanel.setAvailableObjects(objects);
            }
        };
        StsToolkit.runRunnable(runnable);
    }

    public void runGetProjectObjects(StsOSWellObjectTransferPanel _transferPanel, StsProgressPanel _progressPanel)
    {
        final StsOSWellObjectTransferPanel transferPanel = _transferPanel;
        final StsProgressPanel progressPanel = _progressPanel;
        System.out.println("Before Thread:Datastore = " + datastoreName + " Project = " + projectName);
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                Object[] objects = getProjectObjects(progressPanel, null);
                transferPanel.setAvailableObjects(objects);
            }
        };
        StsToolkit.runRunnable(runnable);
    }

    public Object[] getProjectObjects(StsProgressPanel progressPanel, String[] wellboreList)
    {
        System.out.println("After Thread:Datastore = " + datastoreName + " Project = " + projectName);
    	if (projectName == null || datastoreName == null || projectName.length() == 0)
    		return null;

    	// construct a datastoreName - projectName string to set the query scope
    	String[] dsProj = new String[1];
    	dsProj[0] = this.datastoreName + " - " + this.projectName;

        System.out.println("getProjectObjects:Setting query scope to:" + dsProj[0]);
    	ospImport.setQueryScope(dsProj);

        StsOSWell[] wells = null;

        // Determine how many wells are in project first -- OpenSpirit does not have a way to query the number of wells.
        //System.out.println("Calling ospImport.getNumberOfWells - Time:" + System.currentTimeMillis());
        //int numWells = ospImport.getNumberOfWells(projectName);

        if(wellboreList == null)
            // Get all wells in project
            wells = ospImport.getAllProjectWells(projectName, projectLimit, progressPanel);
        else
            // Get only the specified wells
            wells = ospImport.getWellsFromWellboreLists(wellboreList);
            //wells = ospImport.getProjectWellsInList(projectName, projectLimit, progressPanel, wellboreList);

        // temporary restriction to a wellborelist so we don't get too many wells
//        String[] wbList = new String[1];
//        wbList[0] = "AC025";	// 15 wells but no picks
 //       wbList[0] = "BA022A";	// 34 wells some wells have picks
 //       wells = ospImport.getWellsFromWellboreLists(wbList);

        return wells;
    }


}
