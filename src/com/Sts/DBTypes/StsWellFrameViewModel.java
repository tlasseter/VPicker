package com.Sts.DBTypes;

import com.Sts.MVC.ViewWell.*;
import com.Sts.MVC.*;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Dec 17, 2009
  * Time: 8:38:27 AM
  * To change this template use File | Settings | File Templates.
  */
 public class StsWellFrameViewModel extends StsWellViewModel
 {
     /** well window frame size */
     public Dimension frameSize = null;
     /** well window frame location */
     public Point frameLocation = null;
     /** frame for the view of this well model */
     transient public StsWellWindowFrame wellWindowFrame = null;

     public StsWellFrameViewModel()
     {
     }

     public StsWellFrameViewModel(StsWell well)
     {
         super(well, well.currentModel);
         buildFrameView();
     }

     public void restoreWellView()
     {
         well.wellViewModel = this;
         if(isVisible) buildFrameView();
     }

     public void buildFrameView()
     {
         isVisible = true;
         // save restored values temporarily, as they will get stomped upon during resurrection
         Point oldFrameLocation = null;
         Dimension oldFrameSize = null;
         int zoomOldLevel = zoomCurrentLevel;
         // double mdepthMin = mdepthMax - frameSize.getHeight()*zoomLevel.pixelsPerUnit;
         int wX = windowX;
         int wY = windowY;
         if (frameLocation != null) oldFrameLocation = new Point(frameLocation);
         if (frameSize != null) oldFrameSize = new Dimension(frameSize);

         layoutWellWindow();

         if (frameLocation != null)
         {
             wX = frameLocation.x;
             wY = frameLocation.y;
         }

         wellWindowFrame = new StsWellWindowFrame(this, new Point(wX, wY));
         actionManager = wellWindowFrame.actionManager;
         initializeWellWindowPanel();
         wellWindowFrame.add(wellWindowPanel);
         if (oldFrameLocation != null) wellWindowFrame.setLocation(oldFrameLocation);
         if (oldFrameSize != null) wellWindowFrame.setSize(oldFrameSize);
         if (oldFrameLocation != null) wellWindowFrame.setLocation(oldFrameLocation);
         restoreZoom(zoomOldLevel);
         well.setWellFrameViewModel(this);
         start();
     }

     public Frame getParentFrame() { return wellWindowFrame; }
     protected Component getCenterComponent() { return wellWindowFrame; }

     public void setFrameSize(Dimension d)
     {
         frameSize = new Dimension(d);
     }

     public void setFrameLocation(Point p)
     {
         frameLocation = new Point(p);
     }

     public void rebuild()
     {
         wellWindowFrame.rebuild();
     }

     public void printWindow()
     {
         wellWindowFrame.printWindow();
     }// this well might be reused after closing wellWindow, so reinitialize wellMin & wellMax

     // but retain color and name.  Remove all isVisible logs and tracks
           public boolean closeWindow()
           {
               isVisible = false;
               if (wellWindowFrame != null)
               {
                   wellWindowFrame.dispose();
               }
               wellWindowFrame = null;
               well.close();
               for (StsWellView wellView : wellViews)
                   wellView.glPanel.destroy();
               return true;
           }
  
           public void start()
           {
               wellWindowFrame.start();
           }

           public StsWellWindowFrame getWellWindowFrame()
           {
               return wellWindowFrame;
           }

           public void setWellWindowFrame(StsWellWindowFrame frame)
           {
               wellWindowFrame = frame;
           }
       }
