package com.Sts.Actions.Wizards.MakeMovie;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsMakeMovieDefinePanel extends StsJPanel implements StsFileTransferObjectFaceNew
{
    private StsMakeMovieWizard wizard;
    private StsMakeMovieDefine wizardStep;

    boolean start = false;
    Image image = null;
    
    StsGroupBox captureBox = new StsGroupBox("Capture Movie Frames");
    StsToggleButton startStopBtn = null;

    int frameRate = 2;
    //boolean highResolution = false;
    String movieName = null;
    StsGroupBox outputBox = new StsGroupBox("Define Output");
    StsIntFieldBean frameRateBean = null;
    StsStringFieldBean nameBean = null;
    //StsBooleanFieldBean highResBean = null;
    private StsFileTransferPanel filePanel;

    StsGroupBox imageBox = new StsGroupBox("Image Preview");
    JPanel imagePanel = null;
    JScrollPane scrollPane = null;

    private String outputDirectory = null;
    private StsFile currentSelectedFile = null;

    public StsMakeMovieDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMakeMovieWizard)wizard;
        this.wizardStep = (StsMakeMovieDefine)wizardStep;

        try
        {
            outputDirectory = wizard.getModel().getProject().getRootDirString() + wizard.getModel().getProject().getMediaDirString();
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
    	scrollPane = new JScrollPane();
        startStopBtn = new StsToggleButton("Capture Frames","Capture all frames during visualization.", this, "startStop");
        frameRateBean = new StsIntFieldBean(this,"frameRate",1,20,"Frame Rate(f/sec):",true);
        nameBean = new StsStringFieldBean(this, "movieName", true, "Name:");
        //highResBean = new StsBooleanFieldBean(this, "highRes", "Use High Resolution Images Only", false);
        long time = System.currentTimeMillis();
        String timeS = String.valueOf(time);
        timeS = timeS.substring(0,timeS.length()-3);
        movieName = "S2SMovie" + timeS + ".mov";
        nameBean.setValue(movieName);

        String[] filterStrings = new String[] { "jpg" };
        StsFilenameEndingFilter filenameFilter = new StsFilenameEndingFilter(filterStrings);
        filePanel = new StsFileTransferPanel(outputDirectory, filenameFilter, this, 350, 150);
    }
	
    public void setHighRes(boolean val)
    {
    	new StsMessage(wizard.frame,StsMessage.INFO,"Currently only supports JPG format.");
    	return;
    	
    	/*
    	highResolution = val;
        String[] filterStrings = new String[] { "jpg" };    	
    	if(highResolution)
    		filterStrings = new String[] { "bmp" };

    	filePanel.removeSelectedFiles();
        StsFilenameEndingFilter filenameFilter = new StsFilenameEndingFilter(filterStrings);
        filePanel.setFilenameFilter(filenameFilter);
        filePanel.refreshAvailableList();
        */
    }
    
    //public boolean getHighRes() { return highResolution; }
    public void initialize()
    {

    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0f;        
        captureBox.addEndRow(startStopBtn);
        //captureBox.addEndRow(highResBean);
        add(captureBox);

        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0f;
        imagePanel = new JPanel();
        imagePanel.setBackground(Color.WHITE);        
        scrollPane.getViewport().add(imagePanel);
        addEndRow(scrollPane);

        gbc.weighty = 0.0f;
        gbc.fill = gbc.HORIZONTAL;
        outputBox.addEndRow(nameBean);
        outputBox.addEndRow(frameRateBean);
        gbc.gridwidth = 2;
        addEndRow(filePanel);
        addEndRow(outputBox);
    }
    public int getFrameRate() { return frameRate; }
    public String getMovieName() { return movieName; }
    public void startStop()
    {
        wizard.getModel().win3d.outputMovie();
        filePanel.refreshAvailableList();
    }
    public boolean isCapturing()
    {
    	return wizard.getModel().win3d.captureMovie();
    }

    public void setFrameRate(int fr) { frameRate = fr; }
    public void setMovieName(String name) { movieName = name; }

    public boolean hasDirectorySelection() { return false;  }
    public boolean hasReloadButton() { return false;  }
    public void setReload(boolean reload) {}
    public boolean getReload()  { return true; }
    public boolean hasArchiveItButton() { return false; }
    public void setArchiveIt(boolean archive) { }
    public boolean getArchiveIt() { return false; }
    public boolean hasOverrideButton() { return false; }
    public void setOverrideFilter(boolean override) {}
    public boolean getOverrideFilter() { return false; }

    public void fileSelected(StsAbstractFile selectedFile)
    {
        if(selectedFile == null)
        {
            currentSelectedFile = null;
        }
        if (selectedFile != currentSelectedFile)
        {
            currentSelectedFile = (StsFile)selectedFile;
            showImage(currentSelectedFile);
        }
    }
    
    public void availableFileSelected(StsAbstractFile selectedFile)
    {
        showImage((StsFile)selectedFile);
    }
    
    public void showImage(StsFile imageFile)
    {
    	imagePanel.removeAll();
		ImageIcon image = new ImageIcon(imageFile.getPathname());
		JLabel label = new JLabel( image );
		imagePanel.add(label);
		validate();
    }
	
    public void addFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile((StsFile)files[n], outputDirectory);
        fileSelected(files[0]);
    }

    public void addAllFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.addFile((StsFile)files[n], outputDirectory);
        fileSelected(files[files.length-1]);
    }

    public void removeFiles(StsAbstractFile[] files)
    {
        for (int n = 0; n < files.length; n++)
            wizard.removeFile((StsFile)files[n]);
    }

    public void removeAllFiles()
    {
        wizard.removeFiles();
    }
}