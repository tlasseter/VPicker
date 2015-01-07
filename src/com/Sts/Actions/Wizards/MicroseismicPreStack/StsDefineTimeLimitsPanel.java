package com.Sts.Actions.Wizards.MicroseismicPreStack;

import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineTimeLimitsPanel extends StsJPanel implements StsSelectRowNotifyListener
{
    private StsMicroPreStackWizard wizard;
    private StsDefineTimeLimits wizardStep;
    private boolean useHeader = true;
    private long timeLong = 0L;
    
    StsGroupBox criteriaBox = new StsGroupBox("Define Time Limits");
    StsDateFieldBean startBean;
    StsLongFieldBean durationBean;
	StsBooleanFieldBean useHeadersBean;
	
    protected StsSEGYFormatRec[] allRecords = null;
	protected StsSEGYFormatRec[] requiredRecords;
	
    StsGroupBox coordinateBox = new StsGroupBox("Coordinates");
    StsIntFieldBean yByteLocationBean;
    StsIntFieldBean xByteLocationBean;
    StsBooleanFieldBean coordinatesBean;
    
    StsGroupBox headerBox = new StsGroupBox("Header Values");
    StsDoubleFieldBean yCheckBean;
    StsDoubleFieldBean xCheckBean;
    StsLongFieldBean startCheckBean;
    
    public StsTablePanelNew statusTablePanel = null;
    public StsProgressPanel progressPanel;
    
    public StsDefineTimeLimitsPanel(StsMicroPreStackWizard wizard, StsDefineTimeLimits wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildTablePanel();
        constructBeans();
    }
    
    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "xOrigin", "yOrigin", "bornDateString", "deathDateString"};
        String[] columnTitles = {"Name", "Start X", "Start Y", "Start Time", "End Time"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles);
        statusTablePanel = new StsTablePanelNew(tableModel);
        statusTablePanel.setLabel("Files");
        statusTablePanel.setSize(400, 100);
        statusTablePanel.initialize();
    }
    
	private void constructBeans()
	{
        startBean = new StsDateFieldBean(wizard, "startTime", true, "Start Time/Date:");
        startBean.setToolTipText("Specify the minimum start date of all selected gathers.");
        durationBean = new StsLongFieldBean(wizard, "gatherDuration", true, "Duration(ms):");
        durationBean.setToolTipText("Specify the length of a gather in milliseconds.");        
        useHeadersBean = new StsBooleanFieldBean(wizard, "useHeaders", "Read Start Time from Headers", false);
        useHeadersBean.setToolTipText("Look up the start times in the SegY trace headers.");
        
        coordinatesBean = new StsBooleanFieldBean(wizard, "coordinates", "Read Coordinates from Headers", false);
        coordinatesBean.setToolTipText("Look up the coordinates from the trace headers.");        
        xByteLocationBean = new StsIntFieldBean(wizard, "xByteLocation", 0, 240, "Receiver X Coordinate Header Location:", true);
        xByteLocationBean.setToolTipText("Specify the byte position in the trace header to read the receivers X coordiante.");       
        yByteLocationBean = new StsIntFieldBean(wizard, "yByteLocation", 0, 240, "Receiver Y Coordinate Header Location:", true);
        yByteLocationBean.setToolTipText("Specify the byte position in the trace header to read the receivers Y coordiante.");
        progressPanel = StsProgressPanel.constructor(5, 50);
    }
	
    public void buildPanel()
    {
        removeAll();

        gbc.fill = GridBagConstraints.BOTH;
		statusTablePanel.addSelectRowNotifyListener(this);
        add(statusTablePanel);
        
        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.addEndRow(useHeadersBean);
        criteriaBox.add(startBean);
        criteriaBox.add(durationBean);

        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        
        coordinateBox.gbc.fill = gbc.HORIZONTAL;
        coordinateBox.add(coordinatesBean);
        coordinateBox.add(xByteLocationBean);
        coordinateBox.add(yByteLocationBean);
        add(coordinateBox); 
        
        gbc.anchor = GridBagConstraints.SOUTH;
        add(progressPanel);

        wizard.rebuild();        
    }

    public void initialize()
    {
    	setCoordinates(true);
    	coordinatesBean.setSelected(true);
    	setUseHeaders(true);
    	useHeadersBean.setSelected(true);
    	
    	xByteLocationBean.setValue(wizard.getSegyFormat().getXLoc());
    	yByteLocationBean.setValue(wizard.getSegyFormat().getYLoc());
    	long duration = (long)(wizard.getSelectedSegyDatasets()[0].getNSamplesPerTrace() * wizard.getSelectedSegyDatasets()[0].getSampleSpacing());
    	durationBean.setValue(duration);
    	wizard.setGatherDuration(duration);
    	computeStartTime(0);
    	
		wizard.setSkipReanalyzeTraces(true);
        statusTablePanel.replaceRows( wizard.getSegyVolumesList());
        
        buildPanel();
		wizard.setSkipReanalyzeTraces(false);
		
		setCheckFields();
    }
    
    public long computeStartTime(int index)
    {
    	byte[] hdr = wizard.getSelectedSegyDatasets()[index].getTraceHeaderBinary(0);
    	StsSEGYFormat fmt = wizard.getSelectedSegyDatasets()[index].getSegyFormat();
    	long year = (int)fmt.getTraceRec("YEAR").getHdrValue(hdr, wizard.getIsLittleEndian()) + 2000 - 1970;
    	long dayOfYear = (int)fmt.getTraceRec("DAY").getHdrValue(hdr, wizard.getIsLittleEndian());
    	long hours = (int)fmt.getTraceRec("HOUR").getHdrValue(hdr, wizard.getIsLittleEndian());
    	long minutes = (int)fmt.getTraceRec("MINUTE").getHdrValue(hdr, wizard.getIsLittleEndian());
    	long seconds = (int)fmt.getTraceRec("SECOND").getHdrValue(hdr, wizard.getIsLittleEndian());
    	timeLong = (seconds + (minutes*60l) + (hours*3600l) + (dayOfYear*24l*3600l) + (year * 365l * 24l * 3600l)) * 1000l;
		String timeString = StsDateFieldBean.convertToString(timeLong);
    	if(timeLong <= 0)
    	{
    		new StsMessage(wizard.frame, StsMessage.WARNING, "Unable to read start time from trace headers.\n" +
    				" Must set manually. Please specify start time for all files and duration of each.\n" +
    				" We are assuming the files are same duration and back to back.");
    		setUseHeaders(false);
    		useHeadersBean.setSelected(false);
    		timeString = StsDateFieldBean.convertToString(System.currentTimeMillis());
    	}
    	startBean.setValue(timeString);
    	return timeLong;
    }
    
    
    public void setCoordinates(boolean hdr)
    {
    	if(hdr)
    	{
    		xByteLocationBean.setEditable(true);
    		yByteLocationBean.setEditable(true);
    	}
    	else
    	{
    		xByteLocationBean.setEditable(false);
    		yByteLocationBean.setEditable(false);    		
    	}
    } 
    
    public void setUseHeaders(boolean use)
    {
    	if(use)
    	{
    		startBean.setEditable(false);    		
    		computeStartTime(0);    		
    	}
    	else
    	{
    		startBean.setEditable(true);    		
    	}
    }
    
    public boolean setCheckFields()
    {
    	// Set XY Origin of volumes
    	wizard.analyzeHeaders();
    	wizardStep.analyzeGrid();
    	
    	StsSeismicBoundingBox[] boxes = wizard.getSelectedVolumes();
    	statusTablePanel.setSelectAll(true);
    	for(int i=0; i<boxes.length; i++)
    	{   		
    		byte[] hdr = wizard.getSelectedSegyDatasets()[i].getTraceHeader();
    		StsSEGYFormat fmt = wizard.getSelectedSegyDatasets()[i].getSegyFormat();
    		double origin = (double)fmt.getTraceRec(StsSEGYFormat.SHT_X).getHdrValue(hdr, wizard.getIsLittleEndian());
    		boxes[i].setXOrigin(origin);
    		origin = (double)fmt.getTraceRec(StsSEGYFormat.SHT_Y).getHdrValue(hdr, wizard.getIsLittleEndian());    		
    		boxes[i].setYOrigin(origin);
    		
        	// Set Born/Death Time of Volumes
    		long stime = 0l;
    		long etime = 0l;
    		if(wizard.getUseHeaders())
    		{
    			long fileTime = computeStartTime(i);
    			stime = fileTime;
    			etime = fileTime + wizard.getGatherDuration(); 
    		}
    		else
    		{   			
        		stime = wizard.getStartTimeLong() + (long)(wizard.getGatherDuration() * (long)i);
        		etime = wizard.getStartTimeLong() + (long)(wizard.getGatherDuration() * (long)(i+1));    			
    		}
    		boxes[i].setBornDate(StsDateFieldBean.convertToString(stime));
    		boxes[i].setDeathDate(StsDateFieldBean.convertToString(etime));
    	}
    	
    	// Update the table
        statusTablePanel.replaceRows(wizard.getSegyVolumesList());   	
		return true;
    }
   
	public void rowsSelected( int[] indices) {}
}
