package com.Sts.Actions.Wizards.SensorPartition;

import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.StsButton;
import com.Sts.Utilities.StsMath;

import javax.swing.*;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineRangesPanel extends StsJPanel
{
    private StsSensorPartitionWizard wizard;
    private StsDefineRanges wizardStep;

    StsIntFieldBean numPartitionsBean = new StsIntFieldBean();
    int numPartitions = 1;

    StsGroupBox criteriaBox = new StsGroupBox("Define Partition Ranges");
    StsDateFieldBean[] startBeans = null;
    StsDateFieldBean[] endBeans = null;
    StsColorComboBoxFieldBean[] colorListBeans = null;
    StsBooleanFieldBean enableRealtimeBean = new StsBooleanFieldBean();

    long[] startTimes = new long[0];
    long[] endTimes = new long[0];
    int[] colorIdx = new int[0];
    StsSensor sensor = null;
    boolean enableRealtime = true;

    StsGroupBox analyzeBox = new StsGroupBox();
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
    StsButton saveFilterButton;

    public StsDefineRangesPanel(StsSensorPartitionWizard wizard, StsDefineRanges wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        StsGroupBox box = new StsGroupBox();
        numPartitionsBean.initialize(this, "numPartitions", 1, 20, "# of Partitions:", true);
        box.add(numPartitionsBean);

        enableRealtimeBean = new StsBooleanFieldBean(this, "enableRealtime", "Disable Range Validation", true);
        box.add(enableRealtimeBean);
        
        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.NORTH;
        gbc.weighty = 0.0f;
        addEndRow(box);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Adjust time ranges for stages.");
        StsJPanel msgPanel = new StsJPanel();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        StsJPanel btnPanel = new StsJPanel();
        saveFilterButton = new StsButton("Save Filter", "Save the current time ranges as a filter to be selectively applied to sensors", this, "saveAsFilter");

        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addEndRow(saveFilterButton);

        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);

        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        addEndRow(criteriaBox);

        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0f;
        addEndRow(analyzeBox);
    }

    public void initialize()
    {
    	criteriaBox.removeAll();
        configurePanel();
    }

    public void configurePanel()
    {
        Object[] sensors = (Object[])wizard.getSelectedSensors();
        sensor = (StsSensor)sensors[0];
        long born = sensor.getBornDateLong();
        long death = sensor.getDeathDateLong();

        remove(criteriaBox);
        remove(analyzeBox);
        criteriaBox = new StsGroupBox("Define Partition Ranges");

        int oldNumPartitions = startTimes.length;
        if(oldNumPartitions > numPartitions)
        {
            startTimes = (long[])StsMath.trimArray(startTimes, numPartitions);
            endTimes = (long[])StsMath.trimArray(endTimes, numPartitions);
            colorIdx = (int[])StsMath.trimArray(colorIdx, numPartitions);
        }
        else
        {
            startTimes = (long[])StsMath.arrayGrow(startTimes, numPartitions - oldNumPartitions);
            endTimes = (long[])StsMath.arrayGrow(endTimes, numPartitions - oldNumPartitions);
            colorIdx = (int[])StsMath.arrayGrow(colorIdx, numPartitions - oldNumPartitions);
            for(int i=oldNumPartitions; i<numPartitions; i++)
            {
                startTimes[i] = born;
                endTimes[i] = death;
                colorIdx[i] = ((StsDynamicSensorClass)wizard.model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor")).defaultSensorProperties.getColorIndex(i);
            }
        }

        startBeans = new StsDateFieldBean[numPartitions];
        endBeans = new StsDateFieldBean[numPartitions];
    	colorListBeans = new StsColorComboBoxFieldBean[numPartitions];

        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        add(criteriaBox);

        criteriaBox.gbc.anchor = gbc.NORTH;
        criteriaBox.addToRow(new JLabel("Partition"));
        criteriaBox.addToRow(new JLabel("Start Time"));
        criteriaBox.addEndRow(new JLabel("End Time"));

        for(int i=0; i<numPartitions; i++)
    	{
    		startBeans[i] = new StsDateFieldBean(this,"start"+i, StsDateFieldBean.convertToString(startTimes[i]), true, null);
            endBeans[i] = new StsDateFieldBean(this,"end"+i, StsDateFieldBean.convertToString(endTimes[i]), true, null);
            colorListBeans[i] = new StsColorComboBoxFieldBean(this, "color"+i, null, StsColor.colors32);
            colorListBeans[i].setSelectedIndex(colorIdx[i]);

            criteriaBox.gbc.fill = gbc.NONE;
            criteriaBox.gbc.anchor = gbc.NORTH;
    		criteriaBox.addToRow(colorListBeans[i]);
    		criteriaBox.gbc.fill = gbc.HORIZONTAL;
    		criteriaBox.addToRow(startBeans[i]);
    		criteriaBox.addEndRow(endBeans[i]);
    	}

        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        add(analyzeBox);
    	wizard.rebuild();
    }
    public boolean getEnableRealtime() { return enableRealtime; }
    public void setEnableRealtime(boolean val)
    {
        enableRealtime = val;
        if(!enableRealtime)
            verifyTimeRanges();
    }
    public int getNumPartitions() { return numPartitions; }
    public void setNumPartitions(int num)
    {
        numPartitions = num;
        configurePanel();
        verifyTimeRanges();
        updateDisplay();
    }

    public int[] getColorIndices() { return colorIdx; }
    public long[] getStartTimes() { return startTimes; }
    public long[] getEndTimes() { return endTimes; }

    public void updateDisplay()
    {
        long[] times = sensor.getTimeCurve(0).getTimeVectorLongs();
        int[] passed = new int[times.length];
        sensor.setClustering(true);
        for(int i=0; i<times.length; i++)
        {
            passed[i] = -1;
            for(int j=0; j<numPartitions; j++)
            {
                if((times[i] >= startTimes[j]) && (times[i] <= endTimes[j]))
                {
                    passed[i] = colorIdx[j];
                    break;
                }
            }
        }
        sensor.setClusters(passed);
        wizard.getModel().viewObjectRepaint(this, sensor);
    }
    
    public void checkSetStartTime(int index, String timeStg)
    {
        if(index > (numPartitions-1))
            return;

        startTimes[index] = StsDateFieldBean.convertToLong(timeStg);
        verifyTimeRanges();
        criteriaBox.updateBeans();
        updateDisplay();
    }

    public void checkSetEndTime(int index, String timeStg)
    {
        if(index > (numPartitions-1))
            return;

        endTimes[index] = StsDateFieldBean.convertToLong(timeStg);
        verifyTimeRanges();
        criteriaBox.updateBeans();
        updateDisplay();
    }

    public void verifyTimeRanges()
    {
        // Confirm that the input time is not outside sensor time range
        for(int i=0; i<numPartitions; i++)
        {
            if(!enableRealtime)     // All checks are disabled since we don't know what valid ranges will be
            {
                if(startTimes[i] < sensor.getBornDateLong())
                    startTimes[i] = sensor.getBornDateLong();
                if(endTimes[i] > sensor.getDeathDateLong())
                    endTimes[i] = sensor.getDeathDateLong();
                if(startTimes[i] > sensor.getDeathDateLong())
                    startTimes[i] = sensor.getDeathDateLong() - 1;
                if(endTimes[i] < sensor.getBornDateLong())
                    endTimes[i] = sensor.getBornDateLong() + 1;


                // Confirm that startTime is after previous endTime
                if(i >= 1)
                {
                    if(startTimes[i] < endTimes[i-1])
                        startTimes[i] = endTimes[i-1];
                }

                // Confirm that endTime is after startTime
                if(endTimes[i] < startTimes[i])
                    endTimes[i] = startTimes[i];
            }
        }
    }

    public boolean saveAsFilter()
    {
        StsSensorPartitionFilter filter = null;

        filter = new StsSensorPartitionFilter("PartitionFilter");
        filter.addTimeRanges(startTimes, endTimes, colorIdx);

        // Must update the static beans list of filters
        StsObject[] sensors = ((StsSensorClass)wizard.model.getCreateStsClass("com.Sts.DBTypes.StsSensor")).getSensors();
        for(int i=0; i<sensors.length; i++)
            ((StsSensor)sensors[i]).updateFilters();

        msgBean.setText("Partition filter has been successfully saved.");
        return true;
    }

    public String getStart0() { return StsDateFieldBean.convertToString(startTimes[0]); }
    public void setStart0(String timeStg) { checkSetStartTime(0, timeStg); }
    public String getStart1() { return StsDateFieldBean.convertToString(startTimes[1]); }
    public void setStart1(String timeStg) { checkSetStartTime(1, timeStg); }
    public String getStart2() { return StsDateFieldBean.convertToString(startTimes[2]); }
    public void setStart2(String timeStg) { checkSetStartTime(2, timeStg); }
    public String getStart3() { return StsDateFieldBean.convertToString(startTimes[3]); }
    public void setStart3(String timeStg) { checkSetStartTime(3, timeStg); }
    public String getStart4() { return StsDateFieldBean.convertToString(startTimes[4]); }
    public void setStart4(String timeStg) { checkSetStartTime(4, timeStg); }
    public String getStart5() { return StsDateFieldBean.convertToString(startTimes[5]); }
    public void setStart5(String timeStg) { checkSetStartTime(5, timeStg); }
    public String getStart6() { return StsDateFieldBean.convertToString(startTimes[6]); }
    public void setStart6(String timeStg) { checkSetStartTime(6, timeStg); }
    public String getStart7() { return StsDateFieldBean.convertToString(startTimes[7]); }
    public void setStart7(String timeStg) { checkSetStartTime(7, timeStg); }
    public String getStart8() { return StsDateFieldBean.convertToString(startTimes[8]); }
    public void setStart8(String timeStg) { checkSetStartTime(8, timeStg); }
    public String getStart9() { return StsDateFieldBean.convertToString(startTimes[9]); }
    public void setStart9(String timeStg) { checkSetStartTime(9, timeStg); }
    public String getStart10() { return StsDateFieldBean.convertToString(startTimes[10]); }
    public void setStart10(String timeStg) { checkSetStartTime(10, timeStg); }
    public String getStart11() { return StsDateFieldBean.convertToString(startTimes[11]); }
    public void setStart11(String timeStg) { checkSetStartTime(11, timeStg); }
    public String getStart12() { return StsDateFieldBean.convertToString(startTimes[12]); }
    public void setStart12(String timeStg) { checkSetStartTime(12, timeStg); }
    public String getStart13() { return StsDateFieldBean.convertToString(startTimes[13]); }
    public void setStart13(String timeStg) { checkSetStartTime(13, timeStg); }
    public String getStart14() { return StsDateFieldBean.convertToString(startTimes[14]); }
    public void setStart14(String timeStg) { checkSetStartTime(14, timeStg); }
    public String getStart15() { return StsDateFieldBean.convertToString(startTimes[15]); }
    public void setStart15(String timeStg) { checkSetStartTime(15, timeStg); }
    public String getStart16() { return StsDateFieldBean.convertToString(startTimes[16]); }
    public void setStart16(String timeStg) { checkSetStartTime(16, timeStg); }
    public String getStart17() { return StsDateFieldBean.convertToString(startTimes[17]); }
    public void setStart17(String timeStg) { checkSetStartTime(17, timeStg); }
    public String getStart18() { return StsDateFieldBean.convertToString(startTimes[18]); }
    public void setStart18(String timeStg) { checkSetStartTime(18, timeStg); }
    public String getStart19() { return StsDateFieldBean.convertToString(startTimes[19]); }
    public void setStart19(String timeStg) { checkSetStartTime(19, timeStg); }

    public String getEnd0() { return StsDateFieldBean.convertToString(endTimes[0]); }
    public void setEnd0(String timeStg) { checkSetEndTime(0, timeStg); }
    public String getEnd1() { return StsDateFieldBean.convertToString(endTimes[1]); }
    public void setEnd1(String timeStg) { checkSetEndTime(1, timeStg); }
    public String getEnd2() { return StsDateFieldBean.convertToString(endTimes[2]); }
    public void setEnd2(String timeStg) { checkSetEndTime(2, timeStg); }
    public String getEnd3() { return StsDateFieldBean.convertToString(endTimes[3]); }
    public void setEnd3(String timeStg) { checkSetEndTime(3, timeStg); }
    public String getEnd4() { return StsDateFieldBean.convertToString(endTimes[4]); }
    public void setEnd4(String timeStg) { checkSetEndTime(4, timeStg); }
    public String getEnd5() { return StsDateFieldBean.convertToString(endTimes[5]); }
    public void setEnd5(String timeStg) { checkSetEndTime(5, timeStg); }
    public String getEnd6() { return StsDateFieldBean.convertToString(endTimes[6]); }
    public void setEnd6(String timeStg) { checkSetEndTime(6, timeStg); }
    public String getEnd7() { return StsDateFieldBean.convertToString(endTimes[7]); }
    public void setEnd7(String timeStg) { checkSetEndTime(7, timeStg); }
    public String getEnd8() { return StsDateFieldBean.convertToString(endTimes[8]); }
    public void setEnd8(String timeStg) { checkSetEndTime(8, timeStg); }
    public String getEnd9() { return StsDateFieldBean.convertToString(endTimes[9]); }
    public void setEnd9(String timeStg) { checkSetEndTime(9, timeStg); }
    public String getEnd10() { return StsDateFieldBean.convertToString(endTimes[10]); }
    public void setEnd10(String timeStg) { checkSetEndTime(10, timeStg); }
    public String getEnd11() { return StsDateFieldBean.convertToString(endTimes[11]); }
    public void setEnd11(String timeStg) { checkSetEndTime(11, timeStg); }
    public String getEnd12() { return StsDateFieldBean.convertToString(endTimes[12]); }
    public void setEnd12(String timeStg) { checkSetEndTime(12, timeStg); }
    public String getEnd13() { return StsDateFieldBean.convertToString(endTimes[13]); }
    public void setEnd13(String timeStg) { checkSetEndTime(13, timeStg); }
    public String getEnd14() { return StsDateFieldBean.convertToString(endTimes[14]); }
    public void setEnd14(String timeStg) { checkSetEndTime(14, timeStg); }
    public String getEnd15() { return StsDateFieldBean.convertToString(endTimes[15]); }
    public void setEnd15(String timeStg) { checkSetEndTime(15, timeStg); }
    public String getEnd16() { return StsDateFieldBean.convertToString(endTimes[16]); }
    public void setEnd16(String timeStg) { checkSetEndTime(16, timeStg); }
    public String getEnd17() { return StsDateFieldBean.convertToString(endTimes[17]); }
    public void setEnd17(String timeStg) { checkSetEndTime(17, timeStg);}
    public String getEnd18() { return StsDateFieldBean.convertToString(endTimes[18]); }
    public void setEnd18(String timeStg) { checkSetEndTime(18, timeStg);}
    public String getEnd19() { return StsDateFieldBean.convertToString(endTimes[19]); }
    public void setEnd19(String timeStg) { checkSetEndTime(19, timeStg); }

    public StsColor getColor0() { return StsColor.colors32[colorIdx[0]]; }
    public StsColor getColor1() { return StsColor.colors32[colorIdx[1]]; }
    public StsColor getColor2() { return StsColor.colors32[colorIdx[2]]; }
    public StsColor getColor3() { return StsColor.colors32[colorIdx[3]]; }
    public StsColor getColor4() { return StsColor.colors32[colorIdx[4]]; }
    public StsColor getColor5() { return StsColor.colors32[colorIdx[5]]; }
    public StsColor getColor6() { return StsColor.colors32[colorIdx[6]]; }
    public StsColor getColor7() { return StsColor.colors32[colorIdx[7]]; }
    public StsColor getColor8() { return StsColor.colors32[colorIdx[8]]; }
    public StsColor getColor9() { return StsColor.colors32[colorIdx[9]]; }
    public StsColor getColor10() { return StsColor.colors32[colorIdx[10]]; }
    public StsColor getColor11() { return StsColor.colors32[colorIdx[11]]; }
    public StsColor getColor12() { return StsColor.colors32[colorIdx[12]]; }
    public StsColor getColor13() { return StsColor.colors32[colorIdx[13]]; }
    public StsColor getColor14() { return StsColor.colors32[colorIdx[14]]; }
    public StsColor getColor15() { return StsColor.colors32[colorIdx[15]]; }
    public StsColor getColor16() { return StsColor.colors32[colorIdx[16]]; }
    public StsColor getColor17() { return StsColor.colors32[colorIdx[17]]; }
    public StsColor getColor18() { return StsColor.colors32[colorIdx[18]]; }
    public StsColor getColor19() { return StsColor.colors32[colorIdx[19]]; }

     public void setColor0(StsColor color) { colorIdx[0] = StsColor.getColorIndex(color, StsColor.colors32);  updateDisplay(); }
     public void setColor1(StsColor color) { colorIdx[1] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor2(StsColor color) { colorIdx[2] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor3(StsColor color) { colorIdx[3] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor4(StsColor color) { colorIdx[4] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor5(StsColor color) { colorIdx[5] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor6(StsColor color) { colorIdx[6] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor7(StsColor color) { colorIdx[7] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor8(StsColor color) { colorIdx[8] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor9(StsColor color) { colorIdx[9] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor10(StsColor color) { colorIdx[10] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor11(StsColor color) { colorIdx[11] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor12(StsColor color) { colorIdx[12] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor13(StsColor color) { colorIdx[13] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor14(StsColor color) { colorIdx[14] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor15(StsColor color) { colorIdx[15] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor16(StsColor color) { colorIdx[16] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor17(StsColor color) { colorIdx[17] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor18(StsColor color) { colorIdx[18] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }
     public void setColor19(StsColor color) { colorIdx[19] = StsColor.getColorIndex(color, StsColor.colors32);   updateDisplay(); }

}