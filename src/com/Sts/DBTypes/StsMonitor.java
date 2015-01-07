//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsMonitor extends StsMainObject implements StsSelectable, StsTreeObjectI
{
    static protected StsObjectPanel objectPanel = null;
    static public StsFieldBean[] displayFields = null;
    static public final byte TIME = 0;
    static public final byte SIZE = 1;
    static public byte[] pollByTypes  = {TIME, SIZE};
    static public String[] pollByStrings = {"Time", "Size"};

    static public final byte NONE = 0;
    static public final byte FILE = 1;
    static public final byte DIRECTORY = 2;
    static public final byte DB = 3;
    static public final byte DIRECT = 4;
    static public byte[] sourceTypes = {NONE, FILE, DIRECTORY, DB, DIRECT};
    static public String[] sourceStrings = {"None", "File", "Directory", "Database", "Direct Connection" };
    static public String[] filterStrings = new String[]	{"csv", "txt", "Csv", "Txt", "CSV", "TXT"};

    static public final byte SENSOR = 0;
    static public final byte WELL = 1;
    static public final byte LOG = 2;
    byte objectType = SENSOR;

    public StsObject monitorObject = null;
    protected boolean enable = true;
    public String source = null;
    public String filename = null;
    public String source2 = null;
    public String filename2 = null;
    protected byte sourceType = NONE;
    boolean reloadFile = false;

    byte pollBy = TIME;
    boolean replaceEvents = true;
    boolean computeAttributes = false;
    public long lastPollTime = 1L;
    public long lastFileSize = 0;
    public int numPolls = 0;
    public int acceptedChanges = 0;
    public int rejectedChanges = 0;

    transient boolean changed = false;
    transient public String lastPollTimeString = "tst";

    static protected StsBooleanFieldBean enableBean = null;
    static protected StsBooleanFieldBean replaceEventsBean = null;
    static protected StsComboBoxFieldBean pollByBean = null;
    static protected StsBooleanFieldBean reloadFileBean = null;
    static protected StsBooleanFieldBean computeAttributesBean = null;
    static protected StsDateFieldBean lastPollTimeBean = null;
    static protected StsIntFieldBean numPollsBean = null;
    static protected StsIntFieldBean acceptedChangesBean = null;
    static protected StsIntFieldBean rejectedChangesBean = null;

    /** default constructor */
    public StsMonitor()
    {
    }

    public StsMonitor(String name, StsObject object, byte type, String devSource, String devFilename,
                      boolean compute, boolean reload, boolean replace, byte pollBy)
    {
    	super(true, name);
        objectType = type;
    	this.monitorObject = object;
    	this.sourceType = FILE;
        this.source = devSource;
        this.filename = devFilename;
        if(filename != null)
    	   this.source = source + "\\" + filename;
        this.reloadFile = reload;
        this.replaceEvents = replace;
        this.computeAttributes = compute;
        this.pollBy = pollBy;
        this.lastPollTime = 0l;
        setLastPollTimeString();
    	refreshObjectPanel();
    }

    public StsMonitor(String name, StsObject object, String source, byte sourceType, String filename,
                      boolean compute, boolean reload, boolean replace, byte pollBy)
    {
    	super(true, name);
        objectType = SENSOR;
    	this.monitorObject = object;
    	this.sourceType = sourceType;
        this.source = source;
        if(filename != null)
    	   this.source = source + "\\" + filename;
        this.filename = filename;
        this.reloadFile = reload;
        this.replaceEvents = replace;
        this.computeAttributes = compute;
        this.pollBy = pollBy;
        this.lastPollTime = 0l;
        setLastPollTimeString();
    	refreshObjectPanel();
    }

    static public StsMonitor nullMonitorConstructor(String name)
    {
        return new StsMonitor(name, null, null, NONE, null, false, false, true, TIME); //, 0);
    }

    public boolean initialize(StsModel model)
    {
        /*
        if(monitorObject != null)
        {
            if(((StsSensor)monitorObject).getNumValues() > 0)
            {
                lastPollTime = ((StsSensor)monitorObject).getTimeMax();
                acceptedChanges = ((StsSensor)monitorObject).getNumValues();
            }
        }
        */
    	setLastPollTimeString();
    	return true;
    }

    public void setMonitorObject(StsObject object) { this.monitorObject = object; }
    public StsObject getMonitorObject() { return monitorObject; }

    public void setEnable(boolean enable) { this.enable = enable; dbFieldChanged("enable", enable);}
    public boolean getEnable() { return enable; }

    public void setPollBy(byte enable) { this.pollBy = enable; dbFieldChanged("pollBy", pollBy);}
    public byte getPollBy() { return pollBy; }
    public void setPollByString(String pollString) { this.pollBy = getPollByFromString(pollString); dbFieldChanged("pollBy", pollBy);}
    public String getPollByString() { return pollByStrings[pollBy]; }

    public void setSource(String source) { this.source = source; dbFieldChanged("source", source);}
    public String getSource() { return source; }

    public void setSourceTypeString(String typeString) { sourceType = getTypeFromString(typeString); dbFieldChanged("sourceType", sourceType);}
    public String getSourceTypeString() { return sourceStrings[sourceType]; }

    public void setSourceType(byte type) { this.sourceType = type; dbFieldChanged("sourceType", sourceType); }
    public byte getSourceType() { return sourceType; }

    public String toString() { return getName().toString(); }
    public boolean hasChanged() { return changed; }
    public void resetChanged() { changed = false; }
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
        	   enableBean = new StsBooleanFieldBean(this, "enable", "Enable Monitor:");
               pollByBean = new StsComboBoxFieldBean(this, "pollByString", "Poll By:", pollByStrings);
               replaceEventsBean = new StsBooleanFieldBean(this, "replace", "Replace Events", false);
               reloadFileBean = new StsBooleanFieldBean(this, "reload", "Reload File", false);
               computeAttributesBean = new StsBooleanFieldBean(this, "computeAttributes", "Compute Attributes", false);
               computeAttributesBean.setEditable(false);
               lastPollTimeBean = new StsDateFieldBean(StsMonitor.class, "lastPollTime", false, "Last Poll:");
               numPollsBean = new StsIntFieldBean(StsMonitor.class, "numPolls", false, "# Times Polled:");
               acceptedChangesBean = new StsIntFieldBean(StsMonitor.class, "acceptedChanges", false, "# Accepted Events:");
               rejectedChangesBean = new StsIntFieldBean(StsMonitor.class, "rejectedChanges", false, "# Rejected Events:");
               displayFields = new StsFieldBean[]
               {
            		 enableBean, pollByBean, replaceEventsBean, reloadFileBean, computeAttributesBean,
                     new StsStringFieldBean(StsMonitor.class, "source", true, "Source:"),
                     new StsStringFieldBean(StsMonitor.class, "filename", true, "Filename:"),
                     new StsComboBoxFieldBean(StsMonitor.class, "sourceTypeString", "Source Type:", sourceStrings),
                     lastPollTimeBean, numPollsBean, acceptedChangesBean, rejectedChangesBean
               };
           }
           return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields() { return null; }

    public Object[] getChildren() { return new Object[0]; }

    public boolean anyDependencies() { return false; }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void updateProperties() {   }

    public void treeObjectSelected() {  getMonitorClass().selected(this); }

    static public StsMonitorClass getMonitorClass()
    {
        return (StsMonitorClass) currentModel.getCreateStsClass(StsMonitor.class);
    }

    public void setFilename(String name) { this.filename = name; dbFieldChanged("filename", filename);}
    public String getFilename() { return filename; }

    /** remove a monitor from the instance list*/
    public boolean delete()
    {
        return super.delete();
    }

    /** Initialize well even if on uninitialized section. Return true only if
     *  not on section or section is initialized */
    public boolean initialize()
    {
        return true;
    }
    public boolean getComputeAttributes() { return computeAttributes; }
    public void setComputeAttributes(boolean compute)
    {
        computeAttributes = compute;
        dbFieldChanged("computeAttributes", computeAttributes);                  
    }
    public boolean getReload() { return reloadFile; }
    public void setReload(boolean reload)
    {
        reloadFile = reload;
        dbFieldChanged("reloadFile", reloadFile);
    }
    public boolean getReplace() { return replaceEvents; }
    public void setReplace(boolean replace)
    {
        replaceEvents = replace;
        dbFieldChanged("replaceEvents", replaceEvents);
    }
    public boolean canExport() { return false; }
    public boolean export()
    {
        return true;
    }

    static public byte getPollByFromString(String pollByString)
    {
    	for(int i=0; i<pollByStrings.length; i++)
    	{
    		if(pollByStrings[i].equalsIgnoreCase(pollByString))
    			return pollByTypes[i];
    	}
    	return NONE;
    }

    static public byte getTypeFromString(String typeString)
    {
    	for(int i=0; i<sourceStrings.length; i++)
    	{
    		if(sourceStrings[i].equalsIgnoreCase(typeString))
    			return sourceTypes[i];
    	}
    	return NONE;
    }

    public int getNumPolls() { return numPolls; }
    public int getAcceptedChanges() { return acceptedChanges; }
    public int getRejectedChanges() { return rejectedChanges; }

    // Each time the clock is reset when in RT mode, all monitor
    // objects must poll their source for changes.
    public boolean poll()
    {
        int added = 0;
        long currentFileSize = 0;
    	if(!enable || (!currentModel.getProject().isRealtime()) || source == null)
    		return false;


        // Polling may be called multiple times of project time is being set by different sources.
        if(lastPollTime == currentModel.getProject().getProjectTime())
                return false;

    	//System.out.println("Polling for new data.");
    	long time = System.currentTimeMillis();
        numPolls++;

		if(sourceType == DIRECTORY)
        {
            if(!StsToolkit.newDataInDirectory(source, filterStrings, lastPollTime))
            {
                updatePanel(time);
			    return false;
            }
        }
        else
        {
            if(pollBy == TIME)
            {
                if(!StsToolkit.newData(source, lastPollTime))
                {
                    updatePanel(time);
			        return false;
                }
            }
            else
            {
                //System.out.println("Polling by file size....last file size=" + lastFileSize + " file size=" + StsToolkit.getFileSize(source));
                currentFileSize = StsToolkit.getFileSize(source);
                //System.out.println("lastFileSize=" + lastFileSize + " fileSize=" + fileSize);
                if((lastFileSize < currentFileSize))
                    lastFileSize = currentFileSize;
                else if(lastFileSize > currentFileSize)
                {
                    if(!reloadFile)
                    {
                        StsMessageFiles.errorMessage("File size has been reduced and reload file is not selected. Change monitor to reload file.");
                        updatePanel(time);
                        return false;
                    }
                    else
                        lastFileSize = currentFileSize;
                }
                else
                {
                    updatePanel(time);
                    return false;
                }
            }
        }
        long pollTime = System.currentTimeMillis() - time;
        if(pollTime > currentModel.getProject().getTimeUpdateRate())
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Poll rate is faster than access rate. Turning off real-time to allow user reset. Last poll took " + (pollTime/1000) + " seconds. Recommend 2x this speed.");
            currentModel.getProject().stopProjectTime();
            if(pollBy != TIME)
                lastFileSize = currentFileSize;
            return false;
        }
		changed = true;
        //dbFieldChanged("numPolls", numPolls);		// Only persist if new data is ultimately detected.
    	StsMessageFiles.logMessage("New data detected for " + monitorObject.getName() + " at: " + source);

	    // Each object should handle the reading and update of itself
        if(objectType == SENSOR)
            added = ((StsSensor)monitorObject).addNewData(source, sourceType, lastPollTime, computeAttributes, reloadFile, replaceEvents);
        else if((objectType == WELL) || (objectType == LOG))
            added = ((StsLiveWell)monitorObject).addNewData(source, sourceType, lastPollTime, computeAttributes, reloadFile, replaceEvents);            

        if(added > 0)
        {
            if(reloadFile)
                acceptedChanges = added;
            else
        	    acceptedChanges = acceptedChanges + added;
        }
        else
        {
        	rejectedChanges++;
        	dbFieldChanged("rejectedChanges", rejectedChanges);
        }
	    updatePanel(time);
    	return true;
    }

    public String getLastPollTimeString()
    {
    	return lastPollTimeString;
    }
    public void setLastPollTimeString()
    {
    	lastPollTimeString = currentModel.getProject().getTimeDateFormat().format(new Date(lastPollTime));
    	return;
    }
    public String getLastPollTime()
    {
        String timeStg = StsDateFieldBean.convertToString(lastPollTime);
        return timeStg;
    }
    public void setLastPollTime(String polltime) { }

    public void updatePanel(long time)
    {
        lastPollTime = time;
        setLastPollTimeString();
        if(enableBean != null)
        	enableBean.getValueFromPanelObject();
        else
        	return;
        lastPollTimeBean.getValueFromPanelObject();
        numPollsBean.getValueFromPanelObject();
        acceptedChangesBean.getValueFromPanelObject();
        rejectedChangesBean.getValueFromPanelObject();
        reloadFileBean.getValueFromPanelObject();
        replaceEventsBean.getValueFromPanelObject();
        pollByBean.getValueFromPanelObject();

    }
}
