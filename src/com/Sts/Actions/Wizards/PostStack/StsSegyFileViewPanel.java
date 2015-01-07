package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSegyFileViewPanel extends JPanel implements ActionListener, ChangeListener//, ItemListener
{
        private ButtonGroup zGroup = new ButtonGroup();
        private JLabel jLabel1 = new JLabel();
        private StsTablePanel segyTable = new StsTablePanel(); //JTable segyTable = new JTable();
        private JScrollPane jScrollPane1 = new JScrollPane();
        private JLabel jLabel4 = new JLabel();
        private JLabel jLabel5 = new JLabel();
        private JLabel jLabel6 = new JLabel();
        private JComboBox volumeCombo = new JComboBox();
        private GridBagLayout gridBagLayout1 = new GridBagLayout();
        private StsSeismicBoundingBox currentVolume = null;
        private float currentRowNum = -1;
        private StsSliderBean slider = new StsSliderBean();
        private JComboBox headerFormatCombo = new JComboBox();
        private int headerFormat = StsSEGYFormat.INT4;

        Object[] colnames = { "   1", "   5", "   9", "  13", "  17", "  21", "  25", "  29", "  31", "  33",
            "  35", "  37", "  41", "  45", "  49", "  53", "  57", "  61", "  65", "  69", "  71", "  73",
            "  77", "  81", "  85", "  89", "  91", "  93", "  95", "  97", "  99", " 101", " 103", " 105",
            " 107", " 109", " 111", " 113", " 115", " 117", " 119", " 121", " 123", " 125", " 127", " 129",
            " 131", " 133", " 135", " 137", " 139", " 141", " 143", " 145", " 147", " 149", " 151", " 153",
            " 155", " 157", " 159", " 161", " 163", " 165", " 167", " 169", " 171", " 173", " 175", " 177",
            " 179", " 181", " 185", " 189", " 193", " 197", " 201", " 205", " 209", " 213", " 217", " 221",
            " 225", " 229", " 233", " 237" };
        int[] fieldSize = { 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

        Object[] row = new Object[86];

        private DecimalFormat labelFormat = new DecimalFormat("####.#");
        private DecimalFormat intFormat = new DecimalFormat("#");

        private StsWizard wizard;
        private StsWizardStep wizardStep;
        private StsSeismicBoundingBox[] volumes;
        private StsSEGYFormat segyFormat;

        private StsProperties textHeaderFormatTable;
        DefaultTableModel segyTableModel = new DefaultTableModel();

        public StsSegyFileViewPanel(StsWizard wizard, StsWizardStep wizardStep)
        {
            this.wizard = wizard;
            this.wizardStep = wizardStep;
            try
            {
                jbInit();
            }
            catch(Exception ex) {
                ex.printStackTrace();
            }
        }

        void jbInit() throws Exception
        {
            this.setLayout(gridBagLayout1);
            jScrollPane1.setAutoscrolls(true);
            jScrollPane1.setPreferredSize(new Dimension(5000, 5000));

            StsSeismicBoundingBox[] vols = null;
            if(wizard instanceof StsSeismicWizard)
                vols = ((StsSeismicWizard) wizard).getSegyVolumes();
            for(int i=0; i<vols.length; i++)
            {
                volumeCombo.addItem(vols[i].getName());
            }
            volumeCombo.setSelectedIndex(0);
            currentVolume = vols[0];
            currentRowNum = 0;
            volumeCombo.addActionListener(this);

            slider.setIncrementLabel("Step");
            slider.setTextColor(Color.black);
            slider.setValueLabel("Seq Trace:");
//            slider.setPreferredSize(new Dimension(200,40));
            slider.setSelected(true);
            slider.getCheckBoxSlider().setVisible(false);
            initSlider();
            slider.addChangeListener(this);
//            slider.addItemListener(this);

            headerFormatCombo.addItem("IBM Float");
            headerFormatCombo.addItem("IEEE Float");
            headerFormatCombo.addItem("Integer");
            headerFormatCombo.setSelectedIndex(2);
            headerFormatCombo.addActionListener(this);

            segyTable.setPreferredSize(new Dimension(5000, 5000));
            segyTable.setTitle("Traces:");
            segyTable.setFont(new Font("Dialog",3,12));

            jScrollPane1.getViewport().add(segyTable, null);
            this.add(jScrollPane1,  new GridBagConstraints(0, 0, 5, 5, 1.0, 1.0
                ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 8, 0, 17), 0, 0));
            this.add(volumeCombo,  new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
            this.add(headerFormatCombo,  new GridBagConstraints(1, 6, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
            this.add(slider, new GridBagConstraints(2, 6, 3, 1, 1.0, 0.0
                ,GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));

            segyTable.addColumns(colnames);
        }

        public void actionPerformed(ActionEvent e)
        {
            Object source = e.getSource();
            if(source == volumeCombo)
            {
                StsSeismicBoundingBox[] vols = ((StsSeismicWizard) wizard).getSegyVolumes();
                currentVolume = vols[volumeCombo.getSelectedIndex()];
                currentRowNum = 0;
                initSlider();
                setValues();
            }
            else if(source == headerFormatCombo)
            {
//                StsProperties sampleFormats = segyFormat.sampleFormatStrings;
                headerFormat = headerFormatCombo.getSelectedIndex();
                setValues();
            }
        }

        public void initSlider()
        {
            slider.initSliderValues(1, currentVolume.getTotalTraces(), 1, 1);
        }

        public boolean setValues()
        {
            byte[] tbuf = null;
            double val=0.0;

            if(currentVolume == null) return false;

            int startTrace = (int) currentRowNum;
            int nTracesPerLine = currentVolume.getNCols();

            segyTable.removeAllRows();
            segyFormat = currentVolume.getSegyFormat();
            // Segy Trace Table

            if((currentVolume.getNCols() > 0) && (currentVolume.getNCols() < 1000))
                nTracesPerLine = currentVolume.getNCols();
            else
                nTracesPerLine = 1000;

            boolean littleEndian = currentVolume.getIsLittleEndian();
            for (int i = startTrace; i < startTrace + nTracesPerLine; i++)
            {
                tbuf = currentVolume.getTraceHeaderBinary(i);
                int fieldPos = 0;
                for(int j=0; j<row.length; j++)
                {
                    switch(headerFormat)
                    {
                       case 0: // IBM Float
                           val = (float)StsMath.convertIBMFloatBytes(tbuf, fieldPos, littleEndian);
                           row[j] = labelFormat.format(new Float(val));
                           break;
                       case 1: // IEEE Float
                           if(fieldSize[j] == 4)
                               val = (float)Float.intBitsToFloat(StsMath.convertIntBytes(tbuf,fieldPos,littleEndian));
                           else
                               val = (float)Float.intBitsToFloat(StsMath.convertBytesToShort(tbuf,fieldPos,littleEndian));
                           row[j] = labelFormat.format(new Float(val));
                           break;
                       case 2: // Integer
                           if(fieldSize[j] == 4)
                               val = (float)StsMath.convertIntBytes(tbuf, fieldPos,littleEndian);
                           else
                               val = (float)StsMath.convertBytesToShort(tbuf, fieldPos,littleEndian);
                           row[j] = intFormat.format(new Float(val));
                           break;
                        default:
                            ;
                    }
                    fieldPos += fieldSize[j];
                }
                segyTable.addRow(row);
            }
            return true;
        }
/*
        public void itemStateChanged(ItemEvent e)
        {
            currentLineNum = slider.getValue();
            setValues();
            return;
        }
*/
        public void stateChanged(ChangeEvent e)
        {
            StsSliderBean source = (StsSliderBean)e.getSource();
            if(!source.isDraggingSlider())
            {
                currentRowNum = slider.getValue();
                setValues();
            }
            return;
        }

}
