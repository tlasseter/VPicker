package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;
import java.io.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineMonitorPanel extends StsJPanel implements ActionListener
{
    private StsMonitorWizard wizard;
    private StsDefineMonitor wizardStep;
    
    private JFileChooser chooseFile = null; 
    
    byte sourceType = StsMonitor.DIRECTORY;
    String sourceLocation = null;
    String sourceLocation2 = null;
    String filename = null;
    String filename2 = null;
    boolean computeAttributes = false;
    boolean reloadFile = false;
    boolean replaceEvents = true;
    byte pollBy = StsMonitor.TIME;
    String pollByString = null;

    ImageIcon dirIcon = StsIcon.createIcon("dir16x32v2.gif");
    ImageIcon fileIcon = StsIcon.createIcon("file16x32v2.gif");

    StsGroupBox defineBox = new StsGroupBox("Define Monitor");

    StsStringFieldBean sourceLocationBean = null;
    StsStringFieldBean sourceLocation2Bean = null;
    StsComboBoxFieldBean sourceTypeBean = null;
    StsStringFieldBean filenameBean = new StsStringFieldBean();
    StsStringFieldBean filename2Bean = new StsStringFieldBean();
    StsBooleanFieldBean computeAttributesBean = null;
    StsBooleanFieldBean reloadFileBean = null;
    StsBooleanFieldBean replaceOldBean = null;
    StsComboBoxFieldBean pollByBean = null;

    JButton browseButton = new JButton();
    StsGroupBox box = new StsGroupBox();
    JButton browseButton2 = new JButton();
    StsGroupBox box2 = new StsGroupBox();

    
    public StsDefineMonitorPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsDefineMonitor)wizardStep;

        try
        {
            sourceLocation = wizard.getModel().getProject().getRootDirString();
            sourceLocation2 = wizard.getModel().getProject().getRootDirString();
            constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void constructBeans()
    {
        sourceLocationBean = new StsStringFieldBean(this, "sourceLocation", true, "Directory:");
        sourceLocation2Bean = new StsStringFieldBean(this, "sourceLocation2", true, "Directory:");
        sourceTypeBean = new StsComboBoxFieldBean(this, "sourceTypeString", "Source Type:", StsMonitor.sourceStrings);

        pollByBean = new StsComboBoxFieldBean(this, "pollByString", "Poll By:", StsMonitor.pollByStrings);
        pollByBean.setToolTipText("If checked each poll will compare project time to last modified time on file. Alternatively, it will check the file size.");

        computeAttributesBean = new StsBooleanFieldBean(this, "computeAttributes", "Compute Attributes", false);
        computeAttributesBean.setToolTipText("Compute total amplitude and event number attribute.");

        reloadFileBean = new StsBooleanFieldBean(this, "reloadFile", "Reload File", false);
        reloadFileBean.setToolTipText("Reload the entire file on each update.");

        replaceOldBean = new StsBooleanFieldBean(this, "replaceEvents", "Replace Events", false);
        replaceOldBean.setToolTipText("Replace old event with new when time stamp is the same.");

        browseButton.setIcon(dirIcon);
        browseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        browseButton2.setIcon(fileIcon);
        browseButton2.setBorder(BorderFactory.createRaisedBevelBorder());
    }

    public void initialize()
    {
        setSourceTypeString(StsMonitor.sourceStrings[StsMonitor.FILE]);
        sourceTypeBean.setSelectedItem(StsMonitor.sourceStrings[StsMonitor.FILE]);
        if(wizard.getMonitorType() == StsMonitor.WELL)
        {
            sourceTypeBean.setEditable(false);
            filenameBean.initialize(this, "filename", true, "Deviation Filename:");
            filenameBean.setEditable(false);
            filename2Bean.initialize(this, "filename2", true, "Log Filename:");
            filename2Bean.setEditable(false);
            defineBox.add(box2);
        }
        else
        {
            filenameBean.initialize(this, "filename", true, "Sensor Filename:");
            filenameBean.setEditable(false);
            defineBox.remove(box2);
        }
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        defineBox.gbc.fill = gbc.HORIZONTAL;

        box.gbc.fill = gbc.HORIZONTAL;
        box.addEndRow(sourceLocationBean);
        box.addEndRow(filenameBean);
        box.addEndRow(browseButton);

        box2.gbc.fill = gbc.HORIZONTAL;
        box2.addEndRow(sourceLocation2Bean);
        box2.addEndRow(filename2Bean);
        box2.addEndRow(browseButton2);

        defineBox.addEndRow(computeAttributesBean);
        defineBox.addEndRow(reloadFileBean);
        defineBox.addEndRow(replaceOldBean);
        defineBox.addEndRow(pollByBean);

        defineBox.addEndRow(sourceTypeBean);
        defineBox.addEndRow(box);
        defineBox.addEndRow(box2);
        add(defineBox);
        
        browseButton.addActionListener(this);
        browseButton2.addActionListener(this);
    }
    /**
     * Called when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt)
    {
        Object source = evt.getSource();
        File newFile = null;
        try
        {
            if (source == browseButton)
            {
                chooseFile = new JFileChooser(sourceLocation);
                if(sourceType == StsMonitor.FILE)
                    chooseFile.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                else
                    chooseFile.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

                chooseFile.setDialogTitle("Select the data source to monitor.");
                chooseFile.setApproveButtonText("Open");
                while(true)
                {
                    chooseFile.showOpenDialog(null);
                    if(chooseFile.getSelectedFile() == null)
                        break;
                    newFile = chooseFile.getSelectedFile();
                    if(newFile.isDirectory())
                    {
                	    setSourceTypeString(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]);
                        setSourceLocation(newFile.getAbsolutePath());
                        break;
                    }
                    else
                    {
                	    setSourceTypeString(StsMonitor.sourceStrings[StsMonitor.FILE]);
                        setSourceLocation(StsFile.getDirectoryFromPathname(newFile.getPath()));
                        setFilename(newFile.getName());
                        break;
                    }
                } 
        	    updateBeans();
            }
            if (source == browseButton2)
            {
                chooseFile = new JFileChooser(sourceLocation2);
                if(sourceType == StsMonitor.FILE)
                    chooseFile.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                else
                    chooseFile.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

                chooseFile.setDialogTitle("Select the data source to monitor.");
                chooseFile.setApproveButtonText("Open");
                while(true)
                {
                    chooseFile.showOpenDialog(null);
                    if(chooseFile.getSelectedFile() == null)
                        break;
                    newFile = chooseFile.getSelectedFile();

                    if(newFile.isFile())
                    {
                        setSourceLocation2(StsFile.getDirectoryFromPathname(newFile.getPath()));
                        setFilename2(newFile.getName());
                        break;
                    }
                }
        	    updateBeans();
            }
        }
        catch(StsException ex)
        {
            StsException.outputException(ex);
        }
    }
//    public int getFrequency() { return frequency; }
    public String getSourceLocation() { return sourceLocation; }
    public String getSourceLocation2() { return sourceLocation2; }
    public String getSourceTypeString() { return StsMonitor.sourceStrings[sourceType]; }    
    public byte getSourceType() { return sourceType; }

//    public void setFrequency(int fr) { frequency = fr; }
    public void setSourceLocation(String loc) 
    { 
    	sourceLocation = loc; 
    }
    public void setSourceLocation2(String loc)
    {
    	sourceLocation2 = loc;
    }

    public void setSourceTypeString(String typeString) 
    { 
    	if(typeString.equals(StsMonitor.sourceStrings[StsMonitor.DIRECT]) || typeString.equals(StsMonitor.sourceStrings[StsMonitor.DB]))
    	{
    		new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, typeString + " is currently not supported, pick a different type.");
    		setSourceTypeString(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]);
            pollByBean.setEditable(false);
            filenameBean.setEditable(false);
    	}
        else if(typeString.equals(StsMonitor.sourceStrings[StsMonitor.FILE]))
        {
            browseButton.setIcon(fileIcon);
            filenameBean.setEditable(true);
            filename2Bean.setEditable(true);
            pollByBean.setEditable(true);
            wizard.rebuild();
        }
        else
        {
            browseButton.setIcon(dirIcon);
            filenameBean.setEditable(false);
            filename2Bean.setEditable(false);
            pollByBean.setEditable(false);
            setFilename("");
            wizard.rebuild();
        }
    	sourceType = StsMonitor.getTypeFromString(typeString);
    }

    public void setFilename(String name)
    {
        filename = name;
    }
    public String getFilename() { return filename; }
    public void setFilename2(String name)
    {
        filename2 = name;
    }
    public String getFilename2() { return filename2; }

    public void setPollByString(String pollByString)
    {
        pollBy = StsMonitor.getPollByFromString(pollByString);
    }
    public String getPollByString() { return StsMonitor.pollByStrings[pollBy]; }
    public byte getPollBy() { return pollBy; }
    public void setReloadFile(boolean reload)
    {
        reloadFile = reload;
    }
    public boolean getReloadFile() { return reloadFile; }
    public void setReplaceEvents(boolean replace)
    {
        replaceEvents = replace;
    }
    public boolean getReplaceEvents() { return replaceEvents; }
    public void setComputeAttributes(boolean compute)
    {
        computeAttributes = compute;
    }
    public boolean getComputeAttributes() { return computeAttributes; }
    public boolean verifySource(String location)
    {
    	switch(sourceType)
    	{
    		case StsMonitor.FILE:
                if((new File(location).isDirectory()) && filename != null)
    			    return true;
                else
                    return false;
    		case StsMonitor.DIRECTORY:
    			return new File(location).isDirectory();
    		case StsMonitor.DIRECT:        // Not currently supportted
    			return false;
    		case StsMonitor.DB:            // Not currently supportted
    			return false;    			
    	}
    	return false;
    }
}
