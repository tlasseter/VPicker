package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Histogram.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineMathVolumePanel extends JPanel implements ActionListener, ListSelectionListener
{
    private StsWizard wizard;
    private StsDefineMathVolume wizardStep;

    private JPanel selectionPanel = new JPanel();
    private JPanel transferPanel = new JPanel();
    private JButton addBtn = new JButton();
    private JButton addAllBtn = new JButton();
    private JButton removeBtn = new JButton();
    private JButton removeAllBtn = new JButton();

    transient public static final byte ALPHA = 0;
    transient public static final byte AVERAGE = 1;
    transient public static final byte MINIMUM = 2;
    transient public static final byte MAXIMUM = 3;
    transient public static final byte SEGYDATE = 4;
    transient public static String[] SORTBY = {"Alphabetical", "Data Average", "Data Minimum", "Data Maximum", "SegY File Date"};
    JComboBox sortByComboBox = new JComboBox();

//    StsCheckbox pointSetChk = new StsCheckbox();
    private StsSeismicVolume selectedVol = null;
    private int[] selectedIndices;
    private DefaultListModel availableVolsListModel = new DefaultListModel();
    private JScrollPane availableVolsScrollPane = new JScrollPane();

    private DefaultListModel selectedVolsListModel = new DefaultListModel();
    private JScrollPane selectedVolsScrollPane = new JScrollPane();
    private JList selectedVolsList = new JList();

    static final int UP = 0;
    static final int DOWN = 1;

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    private StsModel model = null;
    ButtonGroup modeGroup = new ButtonGroup();
    JButton upBtn = new JButton();
    JButton downBtn = new JButton();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JPanel histPanel = new JPanel();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    StsHistogramPanel histogram = new StsHistogramPanel(StsHistogramPanel.HORIZONTAL);
    JLabel histogramMaxLabel = new JLabel("Maximum");
    JLabel histogramMinLabel = new JLabel("Minimum");
    JLabel histogramMaxPctLabel = new JLabel("Clip %");
    JLabel histogramMinPctLabel = new JLabel("Clip %");
    JLabel maxLabel = new JLabel("Maximum");
    JLabel minLabel = new JLabel("Maximum");
    JLabel histLabel = new JLabel();
    JPanel dataHistogram = new JPanel();

    private int dataCnt[] = new int[255];
    float[] histVals = new float[255];
    int ttlHistogramSamples = 0, minIdx = 0, maxIdx = 254;
    float topClip = 0.0f, btmClip = 0.0f;

    private StsFloatFieldBean dataMinBean;
    private StsFloatFieldBean dataMaxBean;

    private float dataMin = StsParameters.largeFloat;
    private float dataMax = StsParameters.smallFloat;
    private float dataRangeMin = StsParameters.largeFloat;
    private float dataRangeMax = StsParameters.smallFloat;

    JPanel opsPanel = new JPanel();
    private StsButton resetHist = new StsButton();
    JComboBox operatorComboBox = new JComboBox();
    JLabel opLabel = new JLabel("Operator:");

//    private StsCheckbox volumeIdChk = new StsCheckbox();
//    private StsIntFieldBean numSegmentsBean =  null;
//    private int numSegments = 4;

    private DecimalFormat numberFormat = new DecimalFormat("##0.0#");
    JList availableVolsList = new JList();
    GridBagLayout gridBagLayout5 = new GridBagLayout();

    public StsDefineMathVolumePanel(StsWizard wizard, StsDefineMathVolume wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            dataMinBean = new StsFloatFieldBean(this, "dataMin");
            dataMaxBean = new StsFloatFieldBean(this, "dataMax");

            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        if(selectedIndices == null)
        {
            model = wizard.getModel();
            StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);

            int nSeismicVols = volumes.length;
            availableVolsListModel.removeAllElements();
            for(int i=0; i<nSeismicVols; i++)
                availableVolsListModel.addElement(volumes[i].getName());

            volumes = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumes();
            if(volumes != null)
            {
                int nVirtualVols = volumes.length;
                for(int i=0; i<nVirtualVols; i++)
                    availableVolsListModel.addElement(volumes[i].getName());
            }
            //
            // Add operators
            //
            operatorComboBox.removeAllItems();
            for(int i=0; i <= StsMathVirtualVolume.MINIMUM; i++)
                operatorComboBox.addItem(StsMathVirtualVolume.OPERATORS[i]);
            operatorComboBox.setSelectedIndex(StsMathVirtualVolume.MAXIMUM);

            sortByComboBox.removeAllItems();
            for(int i=0; i < SORTBY.length; i++)
                sortByComboBox.addItem(SORTBY[i]);
            sortByComboBox.setSelectedIndex(SEGYDATE);
        }
        resetLimits();
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        selectionPanel.setBorder(BorderFactory.createEtchedBorder());
        selectionPanel.setLayout(gridBagLayout2);
        removeAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeAllBtn.setMaximumSize(new Dimension(100, 20));
        removeAllBtn.setMinimumSize(new Dimension(100, 20));
        removeAllBtn.setPreferredSize(new Dimension(100, 20));
        removeAllBtn.setMargin(new Insets(0, 0, 0, 0));
        removeAllBtn.setText("<< Remove All");
        removeAllBtn.addActionListener(this);
        removeBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        removeBtn.setMaximumSize(new Dimension(100, 20));
        removeBtn.setMinimumSize(new Dimension(100, 20));
        removeBtn.setPreferredSize(new Dimension(100, 20));
        removeBtn.setMargin(new Insets(0, 0, 0, 0));
        removeBtn.setText("< Remove");
        removeBtn.addActionListener(this);
        addAllBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addAllBtn.setMaximumSize(new Dimension(100, 20));
        addAllBtn.setMinimumSize(new Dimension(100, 20));
        addAllBtn.setPreferredSize(new Dimension(100, 20));
        addAllBtn.setMargin(new Insets(0, 0, 0, 0));
        addAllBtn.setText("Add All >>");
        addAllBtn.addActionListener(this);
        transferPanel.setLayout(gridBagLayout1);
        selectedVolsList.setModel(selectedVolsListModel);
        selectedVolsList.addListSelectionListener(this);
        addBtn.setText("Add >");
        addBtn.setFont(new java.awt.Font("Dialog", 0, 11));
        addBtn.setMaximumSize(new Dimension(100, 20));
        addBtn.setMinimumSize(new Dimension(100, 20));
        addBtn.setPreferredSize(new Dimension(100, 20));
        addBtn.setMargin(new Insets(0, 0, 0, 0));
        addBtn.addActionListener(this);
        availableVolsScrollPane.setPreferredSize(new Dimension(150, 110));
        selectedVolsScrollPane.setPreferredSize(new Dimension(150, 110));

        ImageIcon icon = StsIcon.createIcon("upArrow.gif");
        upBtn.setIcon(icon);
        icon = StsIcon.createIcon("downArrow.gif");
        downBtn.setIcon(icon);

        upBtn.addActionListener(this);
        downBtn.addActionListener(this);

        resetHist.setText("Reset Histogram");
        resetHist.setEnabled(true);
        resetHist.setMargin(new Insets(2, 2, 2, 2));
        resetHist.addActionListener(this);
        operatorComboBox.setToolTipText("Select operator to apply");
        opsPanel.setLayout(gridBagLayout5);

        sortByComboBox.setToolTipText("Sort by");
        sortByComboBox.addActionListener(this);

        histPanel.setLayout(gridBagLayout1);
        histLabel.setText("Null Limit Based on Selected Volumes");
        histLabel.setHorizontalAlignment(SwingConstants.CENTER);
        histLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        minLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        minLabel.setText("Minimum:");
        maxLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        maxLabel.setText("Maximum:");
        histogramMinLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinLabel.setForeground(Color.gray);
        histogramMinPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMinPctLabel.setForeground(Color.gray);
        histogramMaxLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxLabel.setForeground(Color.gray);
        histogramMaxPctLabel.setFont(new java.awt.Font("Dialog", 0, 10));
        histogramMaxPctLabel.setForeground(Color.gray);

        JLabel sortByLabel = new JLabel();
        sortByLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        sortByLabel.setText("Sort By:");

        histPanel.setBorder(BorderFactory.createEtchedBorder());
        opsPanel.setBorder(BorderFactory.createEtchedBorder());
        availableVolsList.setModel(availableVolsListModel);
        availableVolsList.addListSelectionListener(this);

        this.add(selectionPanel,      new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 5), 0, 0));

        transferPanel.add(removeBtn,        new GridBagConstraints(0, 3, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addAllBtn,      new GridBagConstraints(0, 2, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));
        transferPanel.add(addBtn,      new GridBagConstraints(0, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 0, 1), 0, 0));
        transferPanel.add(removeAllBtn,         new GridBagConstraints(0, 4, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 4, 0, 1), 0, 0));

        this.add(opsPanel,      new GridBagConstraints(0, 1, 4, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 4, 5, 5), 0, 5));

        opsPanel.add(opLabel,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(4, 5, 0, 0), 0, 0));
        opsPanel.add(operatorComboBox,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 0, 0, 0), 100, 0));
        opsPanel.add(resetHist,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 0, 0, 0), 0, 0));

        this.add(histPanel,  new GridBagConstraints(0, 2, 5, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 4, 5), 0, 0));

        histPanel.add(histLabel, new GridBagConstraints(0, 2, 5, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        histPanel.add(minLabel,  new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0
           ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
        histPanel.add(dataMinBean,  new GridBagConstraints(2, 3, 1, 1, 1.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 2, 0), 0, 0));
        histPanel.add(maxLabel,  new GridBagConstraints(3, 3, 1, 1, 0.0, 0.0
           ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 0, 0), 0, 0));
        histPanel.add(dataMaxBean,  new GridBagConstraints(4, 3, 1, 1, 1.0, 0.0
           ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 0, 0), 0, 0));
        histPanel.add(dataHistogram,  new GridBagConstraints(0, 4, 5, 2, 1.0, 1.0
           ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 0), 0, 0));

        selectionPanel.add(sortByLabel,           new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 4, 0, 1), 0, 0));
        selectionPanel.add(sortByComboBox,           new GridBagConstraints(3, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 4, 0, 1), 0, 0));

        selectionPanel.add(downBtn,   new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 2, 5), 0, 0));
        selectionPanel.add(upBtn,   new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 2, 0), 0, 0));

        selectionPanel.add(availableVolsScrollPane,     new GridBagConstraints(0, 1, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 10, 0), 100, 0));
        availableVolsScrollPane.getViewport().add(availableVolsList, null);
        selectionPanel.add(selectedVolsScrollPane,    new GridBagConstraints(2, 1, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 0, 0, 5), 100, 0));

        selectedVolsScrollPane.getViewport().add(selectedVolsList, null);

        selectionPanel.add(transferPanel,       new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 40, 50));
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        int[] selectedItems = null;

        Object source = e.getSource();

        selectedVolsList.removeListSelectionListener(this);

        // Add selected volumes
        if(source == addBtn)
        {
            selectedItems = availableVolsList.getSelectedIndices();
            for(i=0; i< selectedItems.length; i++)
            {
                if(selectedVolsListModel.indexOf(availableVolsList.getModel().getElementAt(selectedItems[i])) >= 0)
                    continue;
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(selectedItems[i]));
            }
            selectedVolsList.setSelectedIndex(0);
            selectedIndices = selectedVolsList.getSelectedIndices();
            resetLimits();
        }
        // Remove selected volume(s)
        else if (source == removeBtn)
        {
            selectedItems = selectedVolsList.getSelectedIndices();
            for(i= selectedItems.length - 1; i>=0; i--)
                selectedVolsListModel.removeElementAt(selectedItems[i]);
        }
        // Remove all volumes
        else if (source == removeAllBtn)
            selectedVolsListModel.removeAllElements();

        // Add all volumes
        else if (source == addAllBtn)
        {
            selectedVolsListModel.removeAllElements();
            for(i = 0; i < availableVolsList.getModel().getSize(); i++)
                selectedVolsListModel.addElement(availableVolsList.getModel().getElementAt(i));
            selectedVolsList.setSelectedIndex(0);
            selectedIndices = selectedVolsList.getSelectedIndices();
            resetLimits();
        }
        // Up Button
        else if (source == upBtn)
        {
            if(selectedIndices == null)
                new StsMessage(model.win3d, StsMessage.WARNING,"Must select item(s) to move in right list.");
            else
            {
                moveSelectedItem(UP);
            }
        }
        // Down Button
        else if (source == downBtn)
        {
            if(selectedIndices == null)
                new StsMessage(model.win3d, StsMessage.WARNING, "Must select item(s) to move in right list.");
            else
            {
                moveSelectedItem(DOWN);
            }
        }
        else if(source == resetHist)
        {
            dataMin = dataRangeMin;
            dataMax = dataRangeMax;
            resetLimits();
        }
        else if(source == sortByComboBox)
        {
            // Sort the selected volumes as requested.
            sortListBy((byte)sortByComboBox.getSelectedIndex());
        }
        selectedVolsList.addListSelectionListener(this);
    }

    private void moveSelectedItem(int dir)
    {
        StsSeismicVolume[] list = null;
        int i,j, selectedIdx = 0;
        list = getVolumeList();
        if(list == null)
            return;
        switch(dir)
        {
            case UP:
                if(selectedIndices != null)
                {
                    for(j=0; j<selectedIndices.length; j++)
                    {
                        selectedIdx = selectedIndices[j];
                        if(selectedIdx == 0)
                            break;
                        selectedVolsListModel.removeAllElements();
                        for (i = 0; i < selectedIdx - 1; i++)
                            selectedVolsListModel.addElement(list[i].getName());
                        selectedVolsListModel.addElement(list[selectedIdx].getName());
                        selectedVolsListModel.addElement(list[selectedIdx - 1].getName());
                        for (i = selectedIdx + 1; i < list.length; i++)
                            selectedVolsListModel.addElement(list[i].getName());
                    }
                    if (selectedIdx == 0)  break;
                    int[] newSelected = new int[selectedIndices.length];
                    for(j=0; j<selectedIndices.length; j++)
                        newSelected[j] = selectedIndices[j] - 1;
                    selectedVolsList.setSelectedIndices(newSelected);

                    selectedIndices = selectedVolsList.getSelectedIndices();
                }
                break;
            case DOWN:
                if(selectedIndices != null)
                {
                    for(j=selectedIndices.length -1; j>=0; j--)
                    {
                       selectedIdx = selectedIndices[j];
                       if (selectedIdx == (list.length - 1))
                           break;
                       selectedVolsListModel.removeAllElements();
                       for (i = 0; i < selectedIdx; i++)
                           selectedVolsListModel.addElement(list[i].getName());
                       selectedVolsListModel.addElement(list[selectedIdx + 1].getName());
                       selectedVolsListModel.addElement(list[selectedIdx].getName());
                       for (i = selectedIdx + 2; i < list.length; i++)
                           selectedVolsListModel.addElement(list[i].getName());

                    }
                    if (selectedIdx == (list.length - 1))  break;
                    int[] newSelected = new int[selectedIndices.length];
                    for(j=0; j<selectedIndices.length; j++)
                        newSelected[j] = selectedIndices[j] + 1;
                    selectedVolsList.setSelectedIndices(newSelected);

                    selectedIndices = selectedVolsList.getSelectedIndices();
                }

                break;
        }
    }

    public void sortListBy(byte order)
    {
        StsSeismicVolume[] list = null;
        list = getVolumeList();
        selectedVolsListModel.removeAllElements();
        SortData sorter = new SortData();
        sorter.setType(order);
        list = dataSort(list, sorter);

        for(int i=0; i<list.length; i++)
            selectedVolsListModel.addElement(list[i].getName());

        selectedIndices = null;
    }

    private StsSeismicVolume[] dataSort(StsSeismicVolume[] list, Comparator compare)
    {
        Arrays.sort(list, compare);
        return list;
    }

    class SortData implements Comparator
    {
        byte type = SEGYDATE;

        public SortData() {}

        public void setType(byte type)
        {
            this.type = type;
        }

        public int compareTo(Object obj)
        {
            return 0;
        }
        public int compare(Object obj1, Object obj2)
        {
            switch(type)
            {
                case SEGYDATE:
                    return ((StsSeismicVolume)obj1).getName().compareTo(((StsSeismicVolume)obj2).getName());
                case ALPHA:
                    return ((StsSeismicVolume)obj1).getName().compareTo(((StsSeismicVolume)obj2).getName());
                case AVERAGE:
                    if(((StsSeismicVolume)obj1).getDataAvg() < ((StsSeismicVolume)obj2).getDataAvg())
                        return -1;
                    if(((StsSeismicVolume)obj1).getDataAvg() > ((StsSeismicVolume)obj2).getDataAvg())
                        return 1;
                    return 0;
                case MINIMUM:
                    if(((StsSeismicVolume)obj1).getDataMin() < ((StsSeismicVolume)obj2).getDataMin())
                        return -1;
                    if(((StsSeismicVolume)obj1).getDataMin() > ((StsSeismicVolume)obj2).getDataMin())
                        return 1;
                    return 0;
                case MAXIMUM:
                    if(((StsSeismicVolume)obj1).getDataMax() < ((StsSeismicVolume)obj2).getDataMax())
                        return 1;
                    if(((StsSeismicVolume)obj1).getDataMax() > ((StsSeismicVolume)obj2).getDataMax())
                        return -1;
                    return 0;
                default:
                    return 0;
            }
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {
        Object source = e.getSource();
        String volName = null;

        if(source == selectedVolsList)
        {
            selectedIndices = selectedVolsList.getSelectedIndices();
            resetLimits();
        }
        else if(source == availableVolsList)
            selectedIndices = null;

        return;
    }

    public StsSeismicVolume[] getVolumeList()
    {
        String setName = null;
        StsSeismicVolume[] volSet = new StsSeismicVolume[selectedVolsListModel.getSize()];
        for(int i=0; i<selectedVolsListModel.getSize(); i++)
        {
            setName = (String)selectedVolsListModel.getElementAt(i);

            volSet[i] = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, setName);
            if(volSet[i] == null)
            {
                volSet[i] = ((StsVirtualVolumeClass)model.getCreateStsClass(StsVirtualVolume.class)).getVirtualVolumeWithName(setName);
            }
        }
        return volSet;
    }

    private void resetLimits()
    {
        if(selectedIndices == null)
            return;

        clearHistogram();
        dataRangeMin = StsParameters.largeFloat;
        dataRangeMax = StsParameters.smallFloat;
        for(int i=0; i<selectedIndices.length; i++)
        {
            selectedVol = (StsSeismicVolume)model.getObjectWithName(StsSeismicVolume.class, (String) selectedVolsList.getModel().getElementAt(selectedIndices[i]));
            if((dataMin < selectedVol.getDataMin()) || (dataMin == StsParameters.largeFloat))
                dataMin = selectedVol.getDataMin();
            if((dataMax > selectedVol.getDataMax()) || (dataMax == StsParameters.smallFloat))
                dataMax = selectedVol.getDataMax();
            if(dataRangeMin > selectedVol.getDataMin())
                dataRangeMin = selectedVol.getDataMin();
            if(dataRangeMax < selectedVol.getDataMax())
                dataRangeMax = selectedVol.getDataMax();
        }

        int ttlSamples = 0;
        int ttlTopUnusedSamples = 0;
        int ttlBtmUnusedSamples = 0;
        for(int i=0; i<selectedIndices.length; i++)
        {
            selectedVol = (StsSeismicVolume)model.getObjectWithName(StsSeismicVolume.class, (String) selectedVolsList.getModel().getElementAt(selectedIndices[i]));
            ttlSamples = selectedVol.nCols * selectedVol.nRows * selectedVol.nSlices;
            float dataRange = selectedVol.dataMax - selectedVol.dataMin;
            for(int j=0; j<selectedVol.dataHist.length; j++)
            {
                float val = ((float)j/(float)selectedVol.dataHist.length) * dataRange + selectedVol.dataMin;
                int count = (int)((float)selectedVol.dataHist[j]/100.0f * (float)ttlSamples);
                if((val >= dataMin) && (val <= dataMax))
                    accumulateHistogram(val, count);
                else
                {
                    if(val < dataMin)
                        ttlBtmUnusedSamples += count;
                    else
                        ttlTopUnusedSamples += count;
                }
            }
        }
        int total = ttlBtmUnusedSamples + ttlTopUnusedSamples + ttlHistogramSamples;
        topClip = (float)ttlTopUnusedSamples/(float)total * 100.0f;
        btmClip = (float)ttlBtmUnusedSamples/(float)total * 100.0f;

        clearHistogramPanel();
        calculateHistogram();
        displayHistogram(histVals, dataMin, dataMax);

        dataMinBean.setValueAndRange(dataMin, dataRangeMin, dataRangeMax);
        dataMaxBean.setValueAndRange(dataMax, dataRangeMin, dataRangeMax);
    }

    private void clearHistogram()
    {
        for(int i=0; i< 255; i++)
        {
            dataCnt[i] = 0;
            histVals[i] = 0.0f;
        }
        ttlHistogramSamples = 0;
    }

    private void accumulateHistogram(float bindex, int count)
    {
        byte bsamp = 0;
        float scale = 254 / (dataMax - dataMin);
        float scaledFloat = (bindex - dataMin)*scale;
        int scaledInt = Math.round(scaledFloat);
        bsamp = unsignedIntToSignedByte254(scaledInt);

        int index = StsMath.signedByteToUnsignedInt(bsamp);
        if(index > 254) index = 254;
        if(index < 0) index = 0;
        dataCnt[index] = dataCnt[index] + count;
        ttlHistogramSamples += count;
    }

    /** converts an unsigned int to a signedByte value between/including 0 to 254 */
    final public byte unsignedIntToSignedByte254(int i)
    {
        if(i >= 255) i = 254;
        if(i < 0) i = 0;
        return (byte)i;
    }

    private void calculateHistogram()
    {
        for(int i=0; i<255; i++)
            histVals[i] = (float)((float)dataCnt[i]/(float)ttlHistogramSamples)*100.0f;
    }

    public void displayHistogram(float[] histData, float dataMin, float dataMax)
    {
        histogram.updateData(histData, dataMin, dataMax);

        String label = new String(numberFormat.format(btmClip) + "%");
        histogramMinPctLabel.setText(label);
        dataHistogram.add(histogramMinPctLabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        dataHistogram.add(histogram,  new GridBagConstraints(1, 0, 3, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 200, 0));

        label = new String(" " + numberFormat.format(topClip) + "%");
        histogramMaxPctLabel.setText(label);
        dataHistogram.add(histogramMaxPctLabel,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        histogram.repaint();
        dataHistogram.repaint();
        validate();
    }

    public void displayHistogram(StsSeismicVolume volume)
    {
        histogram.updateData(volume.dataHist, volume.getDataMin(), volume.getDataMax());

        String label = new String(numberFormat.format(histogram.getBottomPercentageClipped()) + "%");
        histogramMinPctLabel.setText(label);
        dataHistogram.add(histogramMinPctLabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        dataHistogram.add(histogram,  new GridBagConstraints(1, 0, 3, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 200, 0));

        label = new String(" " + numberFormat.format(histogram.getTopPercentageClipped()) + "%");
        histogramMaxPctLabel.setText(label);
        dataHistogram.add(histogramMaxPctLabel,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 30, 0));

        histogram.repaint();
        validate();
    }

    public void clearHistogramPanel()
    {
        dataHistogram.remove(histogram);
        dataHistogram.remove(histogramMaxLabel);
        dataHistogram.remove(histogramMinLabel);
        dataHistogram.remove(histogramMinPctLabel);
        dataHistogram.remove(histogramMaxPctLabel);
        histogram.clearAll();
    }

    public byte getOperator()
    {
        return (byte)operatorComboBox.getSelectedIndex();
    }

    public float getDataMax() { return dataMax; }
    public void setDataMax(float value)
    {
        dataMax = value;
        if(dataMax > dataRangeMax)
            dataMax = dataRangeMax;
        resetLimits();
    }
    public float getDataMin() { return dataMin; }
    public void setDataMin(float value)
    {
        dataMin = value;
        if(dataMin < dataRangeMin)
            dataMin = dataRangeMin;
        resetLimits();
    }
    public String[] getVolumeNames()
    {
        String[] names = null;

        StsSeismicVolume[] list = getVolumeList();
        names = new String[list.length];
        for(int i=0; i< list.length; i++)
        {
            names[i] = list[i].getName();
        }
        return names;
    }
    public JPanel getHistogramPanel()
    {
        return histPanel;
    }
}
