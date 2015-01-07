package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.util.*;

public class StsAncillaryDataClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassViewSelectable
{
    boolean displayNames = false;
    boolean enableTime = true;
    transient StsSubType[] subTypes = null;

    public static final byte DOCUMENT = 0;
    public static final byte IMAGE = 1;
    public static final byte MULTIMEDIA = 2;
    public static final byte OTHER = 3;
    public static int[] types = { DOCUMENT, IMAGE, MULTIMEDIA, OTHER };

    static String DOC_NODE = "Documents";
    static String IMAGE_NODE = "Images";
    static String MEDIA_NODE = "MultMedia";
    static String OTHER_NODE = "Other";
    public static String[] nodeStrings = { DOC_NODE, IMAGE_NODE, MEDIA_NODE, OTHER_NODE };

    public StsAncillaryDataClass()
    {
        userName = "Docs, Images & Movies";
    }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
           {
           new StsBooleanFieldBean(this, "displayNames", "Names"),
           new StsBooleanFieldBean(this, "enableTime", "Enable Time")
       };

       subTypes = new StsSubType[nodeStrings.length];
       for(int n = 0; n < nodeStrings.length; n++)
           subTypes[n] = new StsSubType(this, nodeStrings[n], (byte)n);
   }

   public void setDisplayNames(boolean displayNames)
   {
       if(this.displayNames == displayNames) return;
       this.displayNames = displayNames;
//       setDisplayField("displayNames", displayNames);
       currentModel.win3dDisplayAll();
   }

   public boolean getDisplayNames() {	return displayNames; }

   public void setEnableTime(boolean enable)
   {
       if(this.enableTime == enable) return;
       this.enableTime = enable;
//       setDisplayField("enableTime", enableTime);
       currentModel.win3dDisplayAll();
   }
   public boolean getEnableTime() {	return enableTime; }

   public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
   {
       if(!isVisible) return;

       int nAncillaryData = getSize();
       if(nAncillaryData == 0) return;
       for(int n = 0; n < nAncillaryData; n++)
       {
           StsAncillaryData data = (StsAncillaryData)getElement(n);
           if((enableTime && data.isAlive(time)) || (!enableTime))
               data.display(glPanel3d, displayNames);
       }
   }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(!isVisible) return;

        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsAncillaryData object = (StsAncillaryData)iter.next();
            object.display(glPanel3d);
        }
    }
    
   public boolean hasDocument(StsWell well)
   {
       int nAncillaryData = getSize();
       if(nAncillaryData == 0)
           return false;
       for(int n = 0; n < nAncillaryData; n++)
       {
           StsAncillaryData doc = (StsAncillaryData)getElement(n);
           if(doc.getWell() == well)
               return true;
       }
       return false;
    }

    public StsSubType[] getSubTypes()
    {
        return subTypes;
    }

    public void createPopupMenu(StsWell well, JMenu menu)
    {
        JMenu subMenu = null;
        StsMenuItem launchItem, exportItem;

        int nAncillaryData = getSize();
        for(int n = 0; n < nAncillaryData; n++)
        {
            StsAncillaryData doc = (StsAncillaryData)getElement(n);
            if(doc.getWell() == well)
            {
                subMenu = new JMenu(doc.getName());
                launchItem = new StsMenuItem();
                launchItem.setMenuActionListener("Launch...", doc, "launch", null);
                subMenu.add(launchItem);
                exportItem = new StsMenuItem();
                exportItem.setMenuActionListener("Export...", doc, "export", null);
                subMenu.add(exportItem);
            }
            if(subMenu != null)
                menu.add(subMenu);
        }
    }
}
