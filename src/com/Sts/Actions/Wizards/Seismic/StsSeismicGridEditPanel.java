package com.Sts.Actions.Wizards.Seismic;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsSeismicGridEditPanel extends StsFieldBeanPanel
{
    public CroppedBoundingBox croppedBoundingBox;
    /** the volumes selected for this pass thru the wizard */
    StsSeismicBoundingBox selectedVolume;

    private boolean rowNumInverted = false;

    private StsGroupBox rangeBox = new StsGroupBox("Volume Crop Box");
    private StsFloatFieldBean rowNumMinBean;
    private StsFloatFieldBean rowNumMaxBean;
    private StsFloatFieldBean colNumMinBean;
    private StsFloatFieldBean colNumMaxBean;
    private StsFloatFieldBean zMinBean;
    private StsFloatFieldBean zMaxBean;
	private JLabel minLabel = new JLabel("Minimum");
	private JLabel maxLabel = new JLabel("Maximum");

    public StsSeismicGridEditPanel()
	{
        try
		{
            croppedBoundingBox = new CroppedBoundingBox();
            constructBeans();
            buildPanel();
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	public void initialize(StsSeismicBoundingBox volume)
	{
        croppedBoundingBox.initialize(volume);
        croppedBoundingBox.initializeBeanValuesAndRanges();
        setPanelObject(croppedBoundingBox);
    }

    private void constructBeans()
    {
        rowNumMinBean = new StsFloatFieldBean(CroppedBoundingBox.class, "rowNumMin", "Inline:");
        rowNumMaxBean = new StsFloatFieldBean(CroppedBoundingBox.class, "rowNumMax");
        colNumMinBean = new StsFloatFieldBean(CroppedBoundingBox.class, "colNumMin", "Crossline:");
        colNumMaxBean = new StsFloatFieldBean(CroppedBoundingBox.class, "colNumMax");
        zMinBean      = new StsFloatFieldBean(CroppedBoundingBox.class, "zMin", "Time/Depth:");
        zMaxBean      = new StsFloatFieldBean(CroppedBoundingBox.class, "zMax");
    }

    private void buildPanel()
	{
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTH;

        rangeBox.gbc.gridx = 1; // leave empty cell
        rangeBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        rangeBox.addToRow(minLabel);
        rangeBox.addEndRow(maxLabel);

        rangeBox.addToRow(rowNumMinBean);
        rangeBox.addEndRow(rowNumMaxBean);

        rangeBox.addToRow(colNumMinBean);
        rangeBox.addEndRow(colNumMaxBean);

        rangeBox.addToRow(zMinBean);
        rangeBox.addEndRow(zMaxBean);

        add(rangeBox);
    }

    public void initializeBeanValuesAndRanges()
    {
        croppedBoundingBox.initializeBeanValuesAndRanges();
    }

    class CroppedBoundingBox extends StsCroppedBoundingBox
    {
        CroppedBoundingBox()
        {
            super(false);
        }

        private void initialize(StsSeismicBoundingBox volume)
        {
            addBoundingBox(volume);

            nSlices = volume.nSlices;
            checkCurrentRangeValues();
            originalBoundingBox = new StsRotatedGridBoundingBox(this);
        }

        /** called when rangeEditPanel is initialized, or one of the range values is changed
         *  in order to check that ranges are ok; if so enable nextButton; otherwise disable it.
         */
        private boolean checkCurrentRangeValues()
        {
            //todo should do a more thorough range check
            return nSlices > 0 && zInc > 0.0f;
        }

        private void initializeBeanValuesAndRanges()
        {
            rowNumInverted = rowNumMin > rowNumMax;
            if(!rowNumInverted)
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMin, rowNumMin, rowNumMax, rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMax, rowNumMin, rowNumMax, rowNumInc);

            }
            else
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMax, rowNumMax, rowNumMin, -rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMin, rowNumMax, rowNumMin, -rowNumInc);
            }

            colNumMinBean.setValueAndRangeFixStep(colNumMin, colNumMin, colNumMax, colNumInc);
            colNumMaxBean.setValueAndRangeFixStep(colNumMax, colNumMin, colNumMax, colNumInc);

            zMinBean.setValueAndRangeFixStep(zMin, zMin, zMax, zInc);
            zMaxBean.setValueAndRangeFixStep(zMax, zMin, zMax, zInc);
        }

        public void setRowNumMin(float rowNumMin)
        {
            if(!rowNumInverted)
                this.rowNumMin = rowNumMin;
            else
                this.rowNumMax = rowNumMin;
        }
        public float getRowNumMin()
        {
            if(!rowNumInverted)
                return rowNumMin;
            else
                return rowNumMax;
        }

        public void setRowNumMax(float rowNumMax)
        {
            if(!rowNumInverted)
                this.rowNumMax = rowNumMax;
            else
                this.rowNumMin = rowNumMax;
        }

        public float getRowNumMax()
        {
            if(!rowNumInverted)
                return rowNumMax;
            else
                return rowNumMin;
        }

        public void setRowNumInc(float value)
        {
//            rowNumMinBean.fixStep(value);
//            rowNumMaxBean.fixStep(value);
        	if(rowNumInc == value)	return;
            rowNumInc = value;
            if(!rowNumInverted)
            {
                rowNumMax = StsMath.floor((originalBoundingBox.rowNumMax - rowNumMin) / rowNumInc) * rowNumInc + rowNumMin;
                rowNumMaxBean.setValue(rowNumMax);
            }
            else
            {
                rowNumMin = StsMath.floor((originalBoundingBox.rowNumMin - rowNumMax) / rowNumInc) * rowNumInc + rowNumMax;
                rowNumMaxBean.setValue(rowNumMin);
            }
        }

        public void setZMin(float value)
        {
            zMin = value;
            // zMax = zMin + (nCroppedSlices-1)*zInc;
            // zMax = Math.min(zMax, originalBoundingBox.zMax);
            // zMaxBean.setValue(zMax);
            // checkCurrentRangeValues();
        }

        public void setZMax(float value)
        {
            zMax = value;
//            zMaxBean.setValue(zMax);
//            zMin = zMax - (nCroppedSlices-1)*zInc;
//            zMinBean.setValue(zMin);
        }
    }

    static public void main(String[] args)
	{
		StsSeismicGridEditPanel panel = new StsSeismicGridEditPanel();
        StsToolkit.createDialog(panel);
    }
}