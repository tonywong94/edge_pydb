#gridplotv2
def gridplotv2(edgetab=None, gallist=None, columnlist=None,
               global_table=None, pixel_table=None, 
               xrange=None, yrange=None, blank=None, plotstyle='image',
               cmap='jet', nx=7, ny=6, dotsize=1, pdfname=None, pct=99,
               allnorm=False, vshow=True, clipedge=False, pad=5, verbose=False, 
               stretch='linear', maxlabel=18, do_cbar=True,
               show_center=False, show_ellipse=False, **kwargs):
    
    from matplotlib.patches import Ellipse

    def get_galaxy_center_pixel(galname, global_table=None, pixel_table=None):
        if global_table is None or pixel_table is None:
            return None
        galdata = global_table[global_table['Name'] == galname]
        if len(galdata) == 0: return None
        try:
            leda_ra = float(galdata['ledaRA'][0])
            leda_dec = float(galdata['ledaDE'][0])
            leda_coord = SkyCoord(ra=leda_ra*u.deg, dec=leda_dec*u.deg, frame='fk5')
        except Exception:
            return None
        galaxy_rows = pixel_table[pixel_table['Name'] == galname]
        if len(galaxy_rows) == 0: return None
        try:
            pixel_coords = SkyCoord(
                ra=np.array(galaxy_rows['ra_abs'])*u.deg,
                dec=np.array(galaxy_rows['dec_abs'])*u.deg,
                frame='fk5'
            )
            separations = leda_coord.separation(pixel_coords)
            nearest_idx = np.argmin(separations)
            ix = int(galaxy_rows['ix'][nearest_idx])
            iy = int(galaxy_rows['iy'][nearest_idx])
            pixel_ra = float(galaxy_rows['ra_abs'][nearest_idx])
            pixel_dec = float(galaxy_rows['dec_abs'][nearest_idx])
            offset_arcsec = separations[nearest_idx].to(u.arcsec).value
        except Exception:
            return None
        return ix, iy, pixel_ra, pixel_dec, offset_arcsec

    def add_ellipse(ax, galname, global_table, pixel_table):
        galdata = global_table[global_table['Name'] == galname]
        galaxy_rows = pixel_table[pixel_table['Name'] == galname]
        if len(galdata) == 0 or len(galaxy_rows) == 0:
            return
        center_info = get_galaxy_center_pixel(galname, global_table, pixel_table)
        if center_info is None:
            return

        ix_center, iy_center, _, _, _ = center_info

        pix_scale = np.abs(galaxy_rows['dec_abs'][1] - galaxy_rows['dec_abs'][0]) * 3600
        Re_arcsec = float(galdata['Re'][0])
        ax_incl_deg = float(galdata['ledaAxIncl'][0])
        pa_deg = float(galdata['ledaPA'][0])

        a_pix = Re_arcsec / pix_scale
        b_pix = a_pix * np.cos(np.radians(ax_incl_deg))
        theta_deg = pa_deg + 90

        # 🔵 ADD DEBUG PRINTS HERE
        print(f"\n=== Ellipse debug for {galname} ===")
        print("ix_center_panel:", ix_center)
        print("iy_center_panel:", iy_center)
        print("a_pix:", a_pix)
        print("b_pix:", b_pix)
        print("theta:", theta_deg)

        ell = Ellipse(
            (ix_center, iy_center),
            width=2*a_pix,
            height=2*b_pix,
            angle=theta_deg,
            edgecolor='red',
            facecolor='none',
            lw=2
        )
        ax.add_patch(ell)
        ax.plot(ix_center, iy_center, marker='o', color='red',
                markersize=14, markeredgewidth=3, markeredgecolor='black')

    match stretch:
        case 'linear': stretch = LinearStretch()
        case 'sqrt': stretch = SqrtStretch()
        case 'log': stretch = LogStretch()
        case _: stretch = LinearStretch()

    if gallist is None and columnlist is None:
        raise TypeError('Either gallist or columnlist must be provided!')
    if isinstance(gallist, str): gallist = [gallist]
    if isinstance(columnlist, str): columnlist = [columnlist]
    if 'norm' in kwargs: kwargs.pop('norm')

    if columnlist is not None and len(columnlist) == 1:
        mode = 'onecol'
        if gallist is None:
            gallist = list(np.unique(edgetab['Name']))
        pagelist = gallist
        if allnorm:
            vmin, vmax = PercentileInterval(pct).get_limits(edgetab[columnlist[0]])
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    elif gallist is not None and len(gallist) == 1:
        mode = 'onegal'
        if columnlist is None:
            columnlist = [c for c in edgetab.colnames if c not in 
                          ['Name','ix','iy','ra_abs','dec_abs','ra_off','dec_off']]
        pagelist = columnlist
        if allnorm:
            vmin, vmax = PercentileInterval(pct).get_limits(
                edgetab[edgetab['Name']==gallist[0]][columnlist])
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    else:
        raise ValueError('Specify either one galaxy or one column to plot')

    pages = int(np.ceil(len(pagelist)/(nx*ny)))
    if pdfname: pp = PdfPages(pdfname)

    for num in range(pages):
        aa, bb = nx*ny*num, nx*ny*(num+1)
        thispage = pagelist[aa:bb]
        fig = plt.figure(figsize=(18,14))
        cbar_done = False

        for i, item in enumerate(thispage):
            ax = plt.subplot(ny,nx,i+1)
            if mode == 'onecol':
                gname = item
                column = columnlist[0]
            else:
                gname = gallist[0]
                column = item
            galtab = (edgetab['Name'] == gname)
            galblank = blank[galtab] if blank is not None else None

            if not np.isnan(edgetab[galtab][column]).all():
                if not allnorm:
                    vmin, vmax = PercentileInterval(pct).get_limits(edgetab[galtab][column])
                    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)

                if plotstyle == 'dot':
                    img, xlims, ylims = dotpatch(edgetab[galtab]['ix'],
                                                 edgetab[galtab]['iy'],
                                                 edgetab[galtab][column],
                                                 blank=galblank, clipedge=clipedge,
                                                 pad=pad, dotsize=dotsize, cmap=cmap,
                                                 norm=norm, axes=ax, **kwargs)
                else:
                    img, xlims, ylims = imarrayplot(edgetab[galtab]['ix'],
                                                   edgetab[galtab]['iy'],
                                                   edgetab[galtab][column],
                                                   blank=galblank, clipedge=clipedge,
                                                   pad=pad, cmap=cmap,
                                                   norm=norm, axes=ax, **kwargs)

                ax.set_xlim(xrange if xrange else xlims)
                ax.set_ylim(yrange if yrange else ylims)

                if show_center:
                    center_info = get_galaxy_center_pixel(gname, global_table, pixel_table)
                    if center_info:
                        ix, iy, _, _, _ = center_info
                        ax.plot(ix, iy, marker='x', color='red', markersize=20, markeredgewidth=2.5)

                if show_ellipse:
                    add_ellipse(ax, gname, global_table, pixel_table)

                if vshow:
                    labelstr = '[{:.3f} .. {:.3f}]'.format(vmin,vmax) if vmax < 1 else '[{:.2f} .. {:.2f}]'.format(vmin,vmax)
                    if hasattr(edgetab[galtab][column], 'unit') and edgetab[galtab][column].unit is not None:
                        if len(edgetab[galtab][column].unit.to_string()) <= maxlabel:
                            labelstr += f" {edgetab[galtab][column].unit:latex_inline}"
                    plt.text(0.04,0.06,labelstr,ha='left',va='center',size='small',transform=ax.transAxes,
                             bbox=dict(boxstyle='square,pad=0.1',facecolor='white',edgecolor='none'))

                if do_cbar and not cbar_done:
                    cax = ax.inset_axes([0.0, -0.07, 1.0, 0.04])
                    cbar = fig.colorbar(img, cax=cax, orientation='horizontal')
                    if allnorm:
                        cbar.ax.tick_params(labelsize='x-small')
                        cbar.ax.tick_params(size=0)
                    else:
                        cbar.set_ticks([])
                    cbar_done = True

            ax.set_aspect('equal')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.text(0.04,0.92, column if mode=='onegal' else gname, ha='left', va='center',
                     transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='none', pad=1))

        fig.subplots_adjust(hspace=0.1, wspace=0.05)
        if pdfname:
            pp.savefig(bbox_inches='tight', pad_inches=0.1)
            plt.close()
        else:
            plt.show()
    if pdfname: pp.close()
