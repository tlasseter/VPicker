package com.Sts.Actions.Wizards.VolumeFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineTypePanel extends StsJPanel
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String KUWAHARA = "Kuwahara";
	public static final String LUM = "LowerUpperMiddle";
	public static final String MSM = "MultiStageMean";
	public static String[] FILTERS = {KUWAHARA};
    private StsVolumeFilterWizard wizard;
    private StsWizardStep wizardStep;

    JPanel jPanel = new JPanel();
    
    StsGroupBox filterGroupBox = new StsGroupBox("Available Filters");
    JCheckBox kuwaharaCheckbox = new JCheckBox();
    JCheckBox LUMCheckbox = new JCheckBox();
   
    StsGroupBox winSizeGroupBox = new StsGroupBox("Window Sizes");
    StsIntFieldBean xSizeBean;
    StsIntFieldBean ySizeBean;
    StsIntFieldBean zSizeBean;
    
    StsGroupBox analysisBox = new StsGroupBox("Define Comparison Criteria");
    ButtonGroup analysisGrp = new ButtonGroup();    
    JRadioButton medianBtn = new JRadioButton("Median");
    JRadioButton meanBtn = new JRadioButton("Mean");



    public StsDefineTypePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsVolumeFilterWizard)wizard;
        this.wizardStep = wizardStep;       

        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
    	meanBtn.setSelected(true);
    }

    void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());
        
        kuwaharaCheckbox.setText("Kuwahara");
        kuwaharaCheckbox.setSelected(true);
//        LUMCheckbox.setText("LowerUpperMiddle");
//        LUMCheckbox.setEnabled(false);
        filterGroupBox.gbc.anchor = GridBagConstraints.WEST;
        filterGroupBox.add(kuwaharaCheckbox);
//        filterGroupBox.add(LUMCheckbox);
        
        analysisGrp.add(medianBtn);
        analysisGrp.add(meanBtn);
        
        analysisBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        analysisBox.gbc.anchor = GridBagConstraints.WEST;        
        analysisBox.addToRow(medianBtn);
        analysisBox.addEndRow(meanBtn);
        
        xSizeBean = new StsIntFieldBean(this, "xSize", 1, 9, "X Size", true);
        xSizeBean.setAlwaysOdd();
        xSizeBean.setValue(3);
        winSizeGroupBox.add(xSizeBean);

        ySizeBean = new StsIntFieldBean(this, "ySize", 1, 9, "Y Size", true);
        ySizeBean.setAlwaysOdd();
        ySizeBean.setValue(3);
        winSizeGroupBox.add(ySizeBean);

        zSizeBean = new StsIntFieldBean(this, "zSize", 1, 3, "Z Size", true);
        zSizeBean.setAlwaysOdd();
        zSizeBean.setValue(1);
        winSizeGroupBox.add(zSizeBean);
        
        this.add(filterGroupBox);
        this.addEndRow(analysisBox);
        this.addEndRow(winSizeGroupBox);

    }

    public ArrayList<String> getAttributes()
    {
        ArrayList<String> list = new ArrayList<String>();
        if (kuwaharaCheckbox.isSelected())
            list.add(StsVolumeFilterConstructor.KUWAHARA);
        if (LUMCheckbox.isSelected())
            list.add(StsVolumeFilterConstructor.LUM);
        return list;
    }
    
    public byte getAnalysisType()
    {
    	if(medianBtn.isSelected())
    		return StsVolumeFilterWizard.ANALYSIS_MEDIAN;
    	else if (meanBtn.isSelected())
    		return StsVolumeFilterWizard.ANALYSIS_MEAN;
    	return StsVolumeFilterWizard.ANALYSIS_NONE;
    }
    
    //public Object getFilterString() { return filterString; }

    //public StsVolumeFilterFace getFilter() { return filter; }

    public void setXSize(int size)
    {
    	xSizeBean.setValue(size);
    }
    public void setYSize(int size)
    {
    	ySizeBean.setValue(size);
    }
    public void setZSize(int size)
    {
    	zSizeBean.setValue(size);
    }

    public int getXSize() { return (int) xSizeBean.getValue(); }
    public int getYSize() { return (int) ySizeBean.getValue(); }
    public int getZSize() { return (int) zSizeBean.getValue(); }

}
