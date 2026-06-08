SUBROUTINE DEALLOCATE_Variables
    USE GlobalVariables
    IMPLICIT NONE
   
#ifdef KE_TURB
if (allocated(vistkg)) DEALLOCATE(vistkg)
if (allocated(visteg)) DEALLOCATE(visteg)
if (allocated(tkg))    DEALLOCATE(tkg)
if (allocated(teg))    DEALLOCATE(teg)

if (allocated(delTkgx)) DEALLOCATE(delTkgx)
if (allocated(delTkgy)) DEALLOCATE(delTkgy)
if (allocated(delTkgz)) DEALLOCATE(delTkgz)
if (allocated(delTegx)) DEALLOCATE(delTegx)
if (allocated(delTegy)) DEALLOCATE(delTegy)
if (allocated(delTegz)) DEALLOCATE(delTegz)
#else
if (allocated(vistkg)) DEALLOCATE(vistkg)
if (allocated(visteg)) DEALLOCATE(visteg)
if (allocated(tkg))    DEALLOCATE(tkg)
if (allocated(teg))    DEALLOCATE(teg)

if (allocated(delTkgx)) DEALLOCATE(delTkgx)
if (allocated(delTkgy)) DEALLOCATE(delTkgy)
if (allocated(delTkgz)) DEALLOCATE(delTkgz)
if (allocated(delTegx)) DEALLOCATE(delTegx)
if (allocated(delTegy)) DEALLOCATE(delTegy)
if (allocated(delTegz)) DEALLOCATE(delTegz)
#endif

#ifdef FIVE_PHASE
if (allocated(uq)) DEALLOCATE(uq)
if (allocated(vq)) DEALLOCATE(vq)
if (allocated(wq)) DEALLOCATE(wq)
if (allocated(phiq)) DEALLOCATE(phiq)
if (allocated(ur)) DEALLOCATE(ur)
if (allocated(vr)) DEALLOCATE(vr)
if (allocated(wr)) DEALLOCATE(wr)
if (allocated(phir)) DEALLOCATE(phir)
if (allocated(us)) DEALLOCATE(us)
if (allocated(vs)) DEALLOCATE(vs)
if (allocated(ws)) DEALLOCATE(ws)
if (allocated(phis)) DEALLOCATE(phis)
#endif
  
  
DEALLOCATE(anu,x,y,z,xcel,ycel,zcel,nod,nc1,nc2,nc3,nc4,vol,dl,sc1x,sc1y,sc1z)
DEALLOCATE(sc2x,sc2y,sc2z,sc3x,sc3y,sc3z,sc4x,sc4y,sc4z,itest,ug,vg,wg,pg,pp )
DEALLOCATE(tg,phig,up,vp,wp,phip,tp,rop,rog,ag,ap,cn,rhs,d24,d2q,cn1,cu,visug)
DEALLOCATE(visvg,viswg,vistg,phi,visup,visvp,viswp,vistp,delrogy,delrogx     )
DEALLOCATE(delphigy,delphigx,delpgz,delvgy,delvgx,deltgy,deltgx,delphigz,delwgz)
DEALLOCATE(delugz,delvgz,delTgz,delwgx,delugx,delrogz,delwgy,delugy,delpgy   )
DEALLOCATE(delupz,delppz,delppx,delppy,delwpx,delwpy,delupx,delupy,delvpy    )
DEALLOCATE(delvpx,deltpy,deltpx,delPgx,delvpz,delTpz,delropy,delropx,delwpz  )
!DEALLOCATE(delropz,cn2,pp_elem,global_pp,neigh_proc) 
DEALLOCATE(neigh_side,dest_ghost,loc_pghost,neigh_glob,cu_send_pos,cu_recv_pos)
!DEALLOCATE(del_send_count,del_recv_count,del_send_pos,del_recv_pos           ) 
!DEALLOCATE(del_send_data,del_recv_data)!,cnd_send_count,cnd_recv_count         )
!DEALLOCATE(cnd_send_data,cnd_recv_data,cnd_send_pos,cnd_recv_pos             )


return
end subroutine DEALLOCATE_Variables

