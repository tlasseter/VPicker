package com.Sts.UI.DataTransfer;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 10, 2008
 * Time: 7:48:32 AM
 * To change this template use File | Settings | File Templates.
 */
/** Interface implemented by datastores which are data sources defined on he StsObjectTransferPanel */
public interface StsDatastoreFace
{
   /** initialize the datastore if necessary */
    public boolean initialize();
    /** sets the name of the selected datstore which has been selected in the StsObjectTransferPanel.browseDatastore() method. */
    public void setDatastoreName(String datastoreName);
    /** returns the names of all available datastores. */
    public String[] getDatastoreNames();
    /** sets new list of available datastore names */
    public void setDatastoreNames(String[] datastoreNames);
     /** returns the name of the datstore which has been selected in the StsObjectTransferPanel.browseDatastore() method. */
    public String getDatastoreName();
    /** Return a complete list of project names available from this datastore. */
    public String[] getProjectNames();
    /** returns name of current project */
    public String getProjectName();
    /** Set the selected projected in the datastore, which should set the available project objects */
    public void setProjectName(String projectName);
    /** Return a list of objects for this datastore/project */
    public Object[] getObjects();
}
