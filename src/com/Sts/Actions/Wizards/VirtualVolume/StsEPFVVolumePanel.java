package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsEPFVVolumePanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;
    /** The volumename is: stemname-filtername
     *  stemname is initialized to volumeName
     *  The user can define a new stemname.
     *  If a new volume is selected, the stemname is changed to the name of this selected volume.
     */
    private String stemname = "";
    /** initialized to stemname-filterName; user can define a stemname.
     *  reset if stemname is changed or if filter is changed.
     *  The user can always override it completely however.
     */
    private String volumeName;

    private StsSeismicVolume volume;

    public static final String KUWAHARA = "Kuwahara";
    //public static final String GAUSSIAN = "Gaussian";
    //public static final String LAPLACIAN = "Laplacian";
    //public static final String MEAN = "Mean";
    //public static final String MEDIAN = "Median";
    //public static final String MAX = "Maximum";
    //public static final String MIN = "Minimum";
    //public static final String VARIANCE = "Variance";

    public static String[] FILTERS = {KUWAHARA};
    
    StsVolumeFilter filter;
    StsGroupBox[] panels;

    StsStringFieldBean volumeNameBean;
    StsComboBoxFieldBean volumeBean;
    StsComboBoxFieldBean filterBean;
    String filterString = KUWAHARA;

    StsJPanel selectPanel;
    StsJPanel filterPanel;
    StsJPanel applyPanel;
    StsButton applyButton;

    StsGroupBox sizePanel;
    StsIntFieldBean xSizeBean;
    StsIntFieldBean ySizeBean;
    StsIntFieldBean zSizeBean;

    public StsEPFVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
    }

    public void setSeismicVolumes(StsSeismicVolume[] seismicVolumes)
    {
        this.seismicVolumes = seismicVolumes;
    }

    public void setVirtualVolumes(StsVirtualVolume[] virtualVolumes)
    {
        this.virtualVolumes = virtualVolumes;
    }

    public void initialize()
    {
        removeAll();

        Object[] volumes = new Object[0];
        volumes = (Object[])StsMath.arrayAddArray(volumes, seismicVolumes);
        volumes = (Object[])StsMath.arrayAddArray(volumes, virtualVolumes);
        volume = (StsSeismicVolume)volumes[0];
        stemname = volume.getName();
        filter = new StsVolumeFilter(filterString);
        volumeName = initializeVolumeName();
        volumeNameBean = new StsStringFieldBean(this, "volumeName","Filter PostStack3d");
        volumeBean = new StsComboBoxFieldBean(this, "volume", "PostStack3d", volumes);
        filterBean = new StsComboBoxFieldBean(this, "filterString","Filter:", FILTERS);

        selectPanel = new StsJPanel();
        //selectPanel.add(volumeNameBean);
        selectPanel.addToRow(volumeBean);
        selectPanel.addEndRow(filterBean);
        add(selectPanel);

        sizePanel = new StsGroupBox("Filter Parameters");
        xSizeBean = new StsIntFieldBean(this, "xSize", 1, 99, "X Size", true);
        xSizeBean.setAlwaysOdd();
        sizePanel.add(xSizeBean);

        ySizeBean = new StsIntFieldBean(this, "ySize", 1, 99, "Y Size", true);
        ySizeBean.setAlwaysOdd();
        sizePanel.add(ySizeBean);

        zSizeBean = new StsIntFieldBean(this, "zSize", 1, 99, "Z Size", true);
        zSizeBean.setAlwaysOdd();
        sizePanel.add(zSizeBean);

        add(sizePanel);

        applyButton = new StsButton("Apply", "Apply filter.", this, "applyFilter");

        rebuildPanel();
    }

    private void rebuildPanel()
    {
        if(filterPanel != null) remove(filterPanel);
        remove(applyButton);

        filter.setName(filterString);
        filterPanel = filter.getFilterPanel();
        if( filterPanel!=null)
            add(filterPanel);

        volumeNameBean.setValue( initializeVolumeName());

        add(applyButton);
        wizard.rebuild();
    }

    public String initializeVolumeName()
    {
        volumeName = stemname +"-"+getFilterString()+filter.getXSize()+filter.getYSize()+filter.getZSize();
        return volumeName;
    }

    public void applyFilter()
    {
        ((StsEPFVVolume)wizardStep).applyFilter(volume, filter, volumeName);
    }

    public String getVolumeName() { return volumeName; }
    public void setVolumeName(String name) { this.volumeName = name;}

    public void setVolume(StsSeismicVolume volume)
    {
        this.volume = volume;
        stemname = volume.getName();
        volumeNameBean.setValue( initializeVolumeName());
    }
    public StsSeismicVolume getVolume() { return volume; }

    public void setFilterString(Object filterString)
    {
        if(this.filterString == filterString) return;
        this.filterString = (String)filterString;
        rebuildPanel();
    }
    public Object getFilterString() { return filterString; }

    public StsVolumeFilterFace getFilter() { return filter; }

    public void setXSize(int size)
    {
        filter.setXSize(size);
        volumeNameBean.setValue( initializeVolumeName());
    }
    public void setYSize(int size)
    {
        filter.setYSize(size);
        volumeNameBean.setValue( initializeVolumeName());
    }
    public void setZSize(int size)
    {
        filter.setZSize(size);
        volumeNameBean.setValue( initializeVolumeName());
    }

    public int getXSize() { return filter.getXSize(); }
    public int getYSize() { return filter.getYSize(); }
    public int getZSize() { return filter.getZSize(); }
}