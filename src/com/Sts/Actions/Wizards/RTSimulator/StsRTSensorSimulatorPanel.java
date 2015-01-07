package com.Sts.Actions.Wizards.RTSimulator;

import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.DateTime.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.text.*;
import java.util.*;

/*
 * ProgressBarDemo.java is a 1.4 application that requires these files:
 *   LongTask.java
 *   SwingWorker.java
 */
public class StsRTSensorSimulatorPanel extends StsJPanel implements ActionListener {
    private String currentDirectory = "C:/Data/Demonstrations/ReservoirMonitoring/Realtime/Events/1G-15_StageF1.csv";
    private String currentFile = "C:/Data/Demonstrations/ReservoirMonitoring/Sensors/1G-15_StageF1.csv";
    private JFileChooser chooseDirectory = null;
    private JFileChooser chooseFile = null;    
    private JFrame frame = null;
    private int delay = 10;
    private String outputTypeString = StsMonitor.sourceStrings[StsMonitor.DIRECTORY];
    private boolean simulateFile = true;
    private boolean enableMultiStage = false;
    private boolean addPastEvent = false;
    private boolean addGarbage = false;
    private int multiStageCriteria = 3500000;  // milliseconds -- less than one hour
    private StsBooleanFieldBean simulateFileBean;
    private StsBooleanFieldBean enableMultiStageBean;
    private StsIntFieldBean delayBean;
    
    protected String ofilename = "1G-15_Stage1.csv";
    protected OutputStream os = null;
    protected BufferedOutputStream bos = null;
    protected ByteArrayOutputStream baos = null;
    protected DataOutputStream dos = null;

    protected InputStream is = null;
    protected BufferedInputStream bis = null;
    protected DataInputStream dis = null;

    private processSensorFile task;

    //JPanel jPanel1 = new JPanel();
    //JPanel jPanel2 = new JPanel();
    StsGroupBox inputBox = new StsGroupBox("Define Input (DD-MM-YY Ordered)");
    StsGroupBox internalBox = new StsGroupBox();
    StsGroupBox outputBox = new StsGroupBox("Define Output");
    StsGroupBox internalBox2 = new StsGroupBox();
    StsGroupBox internalBox3 = new StsGroupBox();
    //JPanel mainPanel = new JPanel();
    JLabel title = new JLabel();
    StsComboBoxFieldBean outputTypeBean = null;    
    JTextField realtimeDirTxt = new JTextField();
    JTextField sensorFileTxt = new JTextField();
    StsDateFieldBean dateBean = new StsDateFieldBean();
    JButton directoryBrowseButton = new JButton();
    JButton fileBrowseButton = new JButton();    
    JButton acceptBtn = new JButton();
    JButton cancelBtn = new JButton();
    JButton pastEventBtn = new JButton();
    JButton garbageBtn = new JButton();
    JButton exitBtn = new JButton();
    SimpleDateFormat format = new SimpleDateFormat("dd/MM/yy kk:mm:ss.SSS");
    SimpleDateFormat outFormat = new SimpleDateFormat("dd/MM/yy,kk:mm:ss.SSS");
    String startDate = new String("01/01/00 00:00:00.000");

	public StsProgressPanel progressPanel = StsProgressPanel.constructorWithCancelButton(5, 50);
	
    public StsRTSensorSimulatorPanel()
    {
        dateBean = new StsDateFieldBean(this, "startDate", "01/01/01 00:00:00.000", true, "Start: ");
        outputTypeBean = new StsComboBoxFieldBean(this, "outputTypeString", "Source Type: ", StsMonitor.sourceStrings);       
        outputTypeBean.setEditable(false);

        task = new processSensorFile(progressPanel);

        delayBean = new StsIntFieldBean(this, "delay", 1, 1000, "Time Multiple:", true);
        delayBean.setValue(10);
        delayBean.setToolTipText("Increase the speed of data flow, one (1) is actual time.");
        simulateFileBean = new StsBooleanFieldBean(this,"simulateFile","Simulate Realtime");
        simulateFileBean.setToolTipText("Replace the time in the file with current time to simulate realtime data flow.");
        enableMultiStageBean = new StsBooleanFieldBean(this,"enableMultiStage","Enable Multiple Stages");
        enableMultiStageBean.setToolTipText("If output to a file, multiple files will be created when time break greater than one hour is detected.");
        this.setLocation(200,200);
        realtimeDirTxt.setText("<supply output directory>");
        realtimeDirTxt.setBackground(Color.WHITE);
        realtimeDirTxt.setToolTipText("Press selection button or type output dirNo or file name followed by return/enter key.");
        sensorFileTxt.setText("<supply sensor file>");
        sensorFileTxt.setBackground(Color.WHITE);   
        sensorFileTxt.setToolTipText("Press selection button or type input file name followed by return/enter key.");
        directoryBrowseButton.setIcon(StsIcon.createIcon("dir16x32v2.gif"));
        directoryBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        fileBrowseButton.setText("Sensor File");
        fileBrowseButton.setBorder(BorderFactory.createRaisedBevelBorder());
        acceptBtn.setText("Accept");
        cancelBtn.setText("Cancel");
        exitBtn.setText("Exit");
        garbageBtn.setText("Add Garbage Line");
        pastEventBtn.setText("Add Event in the Past");

        internalBox.gbc.weightx = 0.0;
        internalBox.gbc.gridwidth = 1;
        internalBox.addToRow(fileBrowseButton);
        internalBox.gbc.weightx = 1.0; 
        internalBox.gbc.gridwidth = 4;
        internalBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        internalBox.addEndRow(sensorFileTxt);

        inputBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        inputBox.gbc.gridwidth = 5;
        inputBox.addEndRow(internalBox);
        inputBox.gbc.gridwidth = 1;
        inputBox.addToRow(simulateFileBean);
        inputBox.addEndRow(enableMultiStageBean);
        inputBox.addToRow(dateBean);

        outputBox.gbc.fill = inputBox.gbc.HORIZONTAL;
        outputBox.addEndRow(outputTypeBean);

        internalBox2.addToRow(directoryBrowseButton);
        internalBox2.gbc.fill = internalBox2.gbc.HORIZONTAL;
        internalBox2.addEndRow(realtimeDirTxt);

        internalBox3.addToRow(garbageBtn);
        internalBox3.addToRow(pastEventBtn);

        outputBox.addEndRow(internalBox2);
        outputBox.addEndRow(delayBean);
        outputBox.addEndRow(internalBox3);

        //add(title);
        gbc.fill = GridBagConstraints.BOTH;        
        addEndRow(inputBox);
        addEndRow(outputBox);
        addEndRow(progressPanel);
        addEndRow(acceptBtn);

        directoryBrowseButton.addActionListener(this);
        fileBrowseButton.addActionListener(this);
        acceptBtn.addActionListener(this);
        garbageBtn.addActionListener(this);
        pastEventBtn.addActionListener(this);
        realtimeDirTxt.addActionListener(this);

        setCurrentDirectory(currentDirectory);
        setCurrentFile(currentFile);
        validateUserInput();
    }

    /**
     * Called when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if (source == fileBrowseButton)
        {
            chooseFile = new JFileChooser(currentFile);
            chooseFile.setCurrentDirectory(new File(currentDirectory));
            chooseFile.setFileSelectionMode(JFileChooser.FILES_ONLY);
            chooseFile.setDialogTitle("Select the input sensor file.");
            chooseFile.setApproveButtonText("Open File");
            while(true)
            {
                chooseFile.showOpenDialog(null);
                if(chooseFile.getSelectedFile() == null)
                    break;
                File newFile = chooseFile.getSelectedFile();
                if (newFile.isFile())
                {
                    setCurrentFile(newFile.getAbsolutePath());
                    break;
                }
                else
                {
                    if(!StsYesNoDialog.questionValue(this,"Must select the file that\n contains the sensor data.\n\n Continue?"))
                        break;
                }
            }
        }        
        else if (source == directoryBrowseButton)
        {
            chooseDirectory = new JFileChooser(currentDirectory);
            chooseDirectory.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            chooseDirectory.setDialogTitle("Select or Enter Desired File or Directory and Press Open");
            chooseDirectory.setApproveButtonText("Open");
            chooseDirectory.showOpenDialog(null);
            File newDirFile = chooseDirectory.getSelectedFile();
            if(newDirFile != null)
            	setCurrentDirectory(newDirFile.getAbsolutePath());
        }
        else if(source == realtimeDirTxt)  // Need carriage return to accept
        { 
        	setCurrentDirectory(realtimeDirTxt.getText());
        }
        else if(source == garbageBtn)
        {
        	addGarbage = true;
        }
        else if(source == pastEventBtn)
        {
        	addPastEvent = true;
        }
        else if(source == acceptBtn)
        {       	
        	if(task.done == true)
        	{
        		if(!validateUserInput())
        			return;
        		StsProgressRunnable progressRunnable = task.getProcessFileRunnable();
        		StsToolkit.runRunnable(progressRunnable);
        	}
        	else
        	{
        		if(StsYesNoDialog.questionValue(frame, "Do you want to re-start?"))
        		{
        			task.progressPanel.cancel();
        			while(!task.done)
        				try {Thread.sleep(delay);} catch(Exception e) { }
        				
            		StsProgressRunnable progressRunnable = task.getProcessFileRunnable();
            		StsToolkit.runRunnable(progressRunnable);
        		}
        	}
        }
    }
    
    public void setOutputTypeString(String typeString) 
    { 
    	if(typeString.equals(StsMonitor.sourceStrings[StsMonitor.DIRECT]) || typeString.equals(StsMonitor.sourceStrings[StsMonitor.DB]))
    	{
    		new StsMessage(this.frame, StsMessage.WARNING, typeString + " is currently not supported, pick a different type.");
    		setOutputTypeString(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]);
    	}
    	outputTypeString = typeString;
    	if(outputTypeString.equalsIgnoreCase(StsMonitor.sourceStrings[StsMonitor.DIRECTORY]))
    		directoryBrowseButton.setIcon(StsIcon.createIcon("dir16x32v2.gif"));
    	else
    		directoryBrowseButton.setIcon(StsIcon.createIcon("file16x32v2.gif"));
    	this.repaint();
    } 
    public String getOutputTypeString() { return outputTypeString; }    
    public boolean validateUserInput()
    {
    	File fileDir = null;
    	try
    	{
            if(!currentFile.equalsIgnoreCase("None"))
            {
    	        fileDir = new File(currentFile);
    		    if(!fileDir.isFile())
    		    {
    			    new StsMessage(this.frame, StsMessage.ERROR, "Selected file is not valid file.");
    			    return false;
    		    }
            }
            
        	fileDir = new File(currentDirectory);
        	if(fileDir.isDirectory()) 
        	{
        		outputTypeBean.doSetValueObject(StsMonitor.sourceStrings[StsMonitor.sourceTypes[StsMonitor.DIRECTORY]]);
        		return true;
        	}
        	else
        	{
        		outputTypeBean.doSetValueObject(StsMonitor.sourceStrings[StsMonitor.sourceTypes[StsMonitor.FILE]]);
        		return true;
        	}        		
    	}
    	catch(Exception e)
    	{
    		new StsMessage(this.frame, StsMessage.ERROR, "Selected input file and/or output directory are not valid.");
    		return false;
    	}
    }
    public boolean getEnableMultiStage() { return enableMultiStage; }
    public void setEnableMultiStage(boolean multi) { this.enableMultiStage = multi; }
    public int getDelay() { return delay; }
    public void setDelay(int delay) { this.delay = delay; }
    public boolean getSimulateFile() { return simulateFile; }
    public void setSimulateFile(boolean time) 
    { 
    	this.simulateFile = time;
    }    
    private void setCurrentDirectory(String directory)
    {
        currentDirectory = directory;
        realtimeDirTxt.setText(currentDirectory);
        validateUserInput();
    }
    public String getStartDate() { return startDate; }
    public void setStartDate(String start) { startDate = start; }
    private void setCurrentFile(String filename)
    {
        currentFile = filename;
        progressPanel.initialize(numberOfLinesInFile(filename));        
        sensorFileTxt.setText(currentFile);
    }
    private int numberOfLinesInFile(String filename)
    {
    	int cnt = 0;
        String line = null;
    	StsFile iFile = StsFile.constructor(filename);
        StsAsciiFile aiFile = new StsAsciiFile(iFile);
        StsSensorFile sensorFile = new StsSensorFile(iFile);
        sensorFile.analyzeFile();
        try
        {
            aiFile.openRead();
            line = aiFile.readLine();  // Header Line.
            line = aiFile.readLine();  // First Data Line.
            while(!sensorFile.getAttributeValues(null,line))
            {
                line = aiFile.readLine();
                if(line == null) return 0;
                continue;
            }
            startDate = format.format(new Date(sensorFile.currentTime));
            if(dateBean != null)
                dateBean.setValue(startDate);
            while(line != null)
            {
        	  cnt++;
              line = aiFile.readLine();
            }
        }
        catch(Exception e)
        {
        	System.out.println("Failed to count number of lines in file:" + filename);
        }
        aiFile.close();
        return cnt - 1;
    }
    
private class processSensorFile {
    private boolean done = true;
    private boolean canceled = false;
    private StsProgressPanel progressPanel;
    
    public processSensorFile(StsProgressPanel progressPanel) {
    	this.progressPanel = progressPanel;
    }
    
    private StsProgressRunnable getProcessFileRunnable()
    {
        StsProgressRunnable progressRunnable = new StsProgressRunnable()
        {
            public void cancel()
            {
                stop();
            }
            public void run()
            {
            	done = false;
            	progressPanel.progressBar.canceled = false;
                readFile();
            }
        };
        return progressRunnable;
    }

    public void stop()
    {
        canceled = true;
        progressPanel.progressBar.canceled = true;
        done = true;
    }

    /**
     * Called from ProgressBarDemo to find out if the task has completed.
     */
    public boolean isDone()
    {
        return done;
    }

    public void readFile()
    {
    	StsFile iFile = null;
        StsAsciiFile aiFile = null;
        String stubname = "Simulated";
        String hline = null;
        String line = null;
        iFile = StsFile.constructor(currentFile);
        aiFile = new StsAsciiFile(iFile);

        try
        {
            long previousTime = ((Date)format.parse(startDate)).getTime() - 1;

            int current = 1;
            StsAbstractFile oFile = null;
            StsAsciiFile outfile = null;
            StsSensorFile sFile = null;

            // Determine the time column(s) from the file
            sFile = new StsSensorFile(iFile);
            sFile.analyzeFile();

            if(iFile.getFilename().contains(".csv"))
                stubname = iFile.getFilename().substring(0,iFile.getFilename().lastIndexOf(".csv"));
            else
                stubname = iFile.getFilename().substring(0,iFile.getFilename().lastIndexOf(".txt"));
            
            aiFile.openRead();
            hline = aiFile.readLine();

            // If simulating realtime, header line is re-ordered.
            if(simulateFile)
            {
                hline = "Date, Time";
                for(int ii=0; ii<sFile.curveNames.length; ii++)
                    hline = hline + "," + sFile.curveNames[ii];
            }
            // If creating file, add header once and then loop through lines.
            String oFilename = currentDirectory;
            if(StsMonitor.getTypeFromString(outputTypeString) == StsMonitor.FILE)
               outfile = createOutputFile(oFilename, hline);

            int stage = 1;
            while(!done)
            {
                // Read or create event
                line = aiFile.readLine();
                sFile.getAttributeValues(null, line, CalendarParser.DD_MM_YY);
                if(line == null)
                   	break;

                // Wait the requested delay.
                if(sFile.currentTime != 0)
                {
                    long waitTime = sFile.currentTime - previousTime;
                    if(waitTime < 0)
                        waitTime = 30000;
                    if(enableMultiStage && (waitTime > multiStageCriteria) && (StsMonitor.getTypeFromString(outputTypeString) == StsMonitor.FILE))
                    {
                        stage++;
                        String directoryName = StsFile.getDirectoryFromPathname(currentDirectory);
                        String fileName = StsFile.getFilenameFromPathname(currentDirectory);
                        String stubName = fileName.substring(0,fileName.indexOf(".csv"));
                        oFilename = directoryName + "/" + stubName + "_Stage" + stage + ".csv";
                        outfile = createOutputFile(oFilename, hline);
                    }
                    long sleepDelay = waitTime/delay;
                    if(previousTime != 0)
                         Thread.sleep(sleepDelay);
                    previousTime = sFile.currentTime;
                }
                // Add a header to every file if directory and open File
              	if(StsMonitor.getTypeFromString(outputTypeString) == StsMonitor.DIRECTORY)
               	{
                   oFilename = currentDirectory + "/" + stubname + "_rt00" + current + ".csv";
                   outfile = createOutputFile(oFilename, hline);
              	}
                outfile.openWriteAppend();

                // Output line
                // --- If simulate realtime, replace time in file with current time.
                long time = sFile.currentTime;
                if(simulateFile)
                {
                    time = System.currentTimeMillis();
                    String timeDate = outFormat.format(new Date(time));
                    line = timeDate;
                    for(int ii=0; ii<sFile.validCurves.length; ii++)
                        line = line + "," + sFile.currentValues[ii];
                }
                //System.out.println("Outputting line to file: " + line);
                outfile.writeLine(line);
               	progressPanel.setValue(current-1);
               	if(StsMonitor.getTypeFromString(outputTypeString) == StsMonitor.DIRECTORY)
                {
                   	progressPanel.appendLine("Output file for time " + format.format(new Date(time)) + ": " + outfile.getFilename());
                }
                else
                   	progressPanel.appendLine("Output line #" + current + " at time " + format.format(new Date(time)) + " to file: " + outfile.getFilename());
                current++;

               	if(progressPanel.isCanceled())
               	{
               		progressPanel.appendLine("Canceling processing... ");
               		done = true;
               	}
                if(addPastEvent)
                {
                    addEventFromPast(outfile, time, sFile);
                }
                if(addGarbage)
                {
                    addGarbageLine(outfile);
                }
                outfile.close();
            }
            if(!canceled)
            {
               	progressPanel.finished();
               	progressPanel.setDescriptionAndLevel("Processing complete", progressPanel.progressBar.INFO);
            }
            else
               	progressPanel.setDescriptionAndLevel("Canceled by user", progressPanel.progressBar.WARNING);
            }
            catch (Exception e) {
            	 progressPanel.setDescriptionAndLevel("Error processing file. Check permissions.", progressPanel.progressBar.ERROR);
                 return;
             }
             finally
             {
                 if(!simulateFile) aiFile.close();
                 done = true;                 
             }
          }
       }

    private void addEventFromPast(StsAsciiFile outfile, long time, StsSensorFile sFile)
    {
        String line = null;
        time = time - 3600000;
        String timeDate = outFormat.format(new Date(time));
        line = timeDate;
        for(int ii=0; ii<sFile.validCurves.length; ii++)
             line = line + "," + sFile.currentValues[ii];
        try { outfile.writeLine(line); }
        catch(Exception ex) { System.out.println("failed to output line in past."); }
        progressPanel.appendLine("Output event from the past - " + timeDate);
        addPastEvent = false;
    }

    private void addGarbageLine(StsAsciiFile outfile)
    {
        String line = "Garbage line...#$ , with , some , trash";
        try { outfile.writeLine(line); }
        catch(Exception ex) { System.out.println("failed to output garbage line."); }
        progressPanel.appendLine("Output garbage line - " + line);
        addGarbage = false;
    }
    
    private StsAsciiFile createOutputFile(String filename, String headerLine)
    {
        StsAsciiFile outfile = null;
        StsAbstractFile oFile = null;

        try
        {
            oFile = StsFile.constructor(filename);
            outfile = new StsAsciiFile(oFile);

            // Truncate file and start on line one.
            outfile.openWrite();
            progressPanel.appendLine("Output header to file: " + outfile.getFilename());
            outfile.writeLine(headerLine);
            outfile.close();
            return outfile;
        }
        catch(Exception ex)
        {
            return null;
        }
    }
}
