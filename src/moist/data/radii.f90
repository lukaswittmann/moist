
module moist_data_radii_legacy
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error
   use mctc_io_convert, only: aatoau
   use mctc_io_symbols, only: symbol_to_number
   use mctc_io_utils, only: to_lower
   implicit none
   private

   integer :: i

   !> Upper bounds for each radii array
   integer, parameter :: max_elem_cpcm = 118
   integer, parameter :: max_elem_smd = 118
   integer, parameter :: max_elem_d3 = 94
   integer, parameter :: max_elem_cosmo = 94
   integer, parameter :: max_elem_bondi = 88
   integer, parameter :: max_elem_rahm = 96
   integer, parameter :: max_elem_gauss = 118

   ! Model tags
   ! TODO: best practice?
   type, public :: radius_type
      integer :: cpcm = 1
      integer :: smd = 2
      integer :: d3 = 3
      integer :: cosmo = 4
      integer :: bondi = 5
      integer :: rahm = 6
      integer :: gauss = 7
   end type radius_type

   type(radius_type), parameter, public :: rad_type = radius_type()

   ! TODO: general handling of different lengths of radii arrays?

   !> CPCM radii (already scaled)
   real(wp), parameter :: cpcm_vdw_rad(max_elem_cpcm) = aatoau*[ &
                          1.300_wp, 1.400_wp*1.17_wp, &
                          1.200_wp*1.17_wp, 0.900_wp*1.17_wp, 1.750_wp*1.17_wp, 2.000_wp, 1.830_wp, 1.720_wp, &
                          1.720_wp, 1.540_wp*1.17_wp, 1.500_wp*1.17_wp, 1.400_wp*1.17_wp, 2.153_wp, 2.200_wp, &
                          1.800_wp*1.17_wp, 2.16_wp, 2.050_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, &
                          [(1.900_wp*1.17_wp, i=1, 10)], &
                          1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, &
                          2.160_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, &
                          [(1.900_wp*1.17_wp, i=1, 10)], &
                          1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 1.900_wp*1.17_wp, 2.320_wp, &
                          1.900_wp*1.17_wp, 2.000_wp*1.17_wp, 2.000_wp*1.17_wp, 2.000_wp*1.17_wp, &
                          [(2.000_wp*1.17_wp, i=1, 14)], [(2.000_wp*1.17_wp, i=1, 6)], &
                          1.720_wp*1.17_wp, 1.720_wp*1.17_wp, 2.000_wp*1.17_wp, &
                          [(2.000_wp*1.17_wp, i=1, 6)], &
                          2.000_wp*1.17_wp, 2.000_wp*1.17_wp, 2.000_wp*1.17_wp, &
                          [(2.000_wp*1.17_wp, i=1, 14)], [(2.000_wp*1.17_wp, i=1, 9)], &
                          [(2.000_wp*1.17_wp, i=1, 6)] &
                          ]

   !> SMD radii
   real(wp), parameter :: smd_vdw_rad(max_elem_smd) = aatoau*[ &
                          1.20_wp, 1.40_wp, &
                          1.82_wp, 1.53_wp, 1.92_wp, 1.85_wp, 1.89_wp, 1.52_wp, 1.73_wp, 1.54_wp, &
                          2.27_wp, 1.73_wp, 1.84_wp, 2.47_wp, 2.12_wp, 2.49_wp, 2.38_wp, 1.88_wp, &
                          2.75_wp, 2.31_wp, &
                          2.16_wp, 1.87_wp, 1.79_wp, 1.89_wp, 1.97_wp, 1.94_wp, 1.92_wp, 1.84_wp, 1.86_wp, 2.10_wp, &
                          1.87_wp, 2.11_wp, 1.85_wp, 1.90_wp, 3.06_wp, 2.02_wp, &
                          3.03_wp, 2.49_wp, &
                          2.19_wp, 1.86_wp, 2.07_wp, 2.09_wp, 2.09_wp, 2.07_wp, 1.95_wp, 2.02_wp, 2.03_wp, 2.30_wp, &
                          1.93_wp, 2.17_wp, 2.06_wp, 2.06_wp, 1.98_wp, 2.16_wp, &
                          3.43_wp, 2.68_wp, &
                          2.40_wp, &
                          2.35_wp, 2.39_wp, 2.29_wp, 2.36_wp, 2.29_wp, 2.33_wp, 2.37_wp, 2.21_wp, 2.29_wp, 2.16_wp, &
                          2.35_wp, 2.27_wp, 2.42_wp, 2.21_wp, &
                          2.12_wp, 2.17_wp, 2.10_wp, 2.17_wp, 2.16_wp, 2.02_wp, 2.09_wp, 2.17_wp, 2.09_wp, &
                          1.96_wp, 2.02_wp, 2.07_wp, 1.97_wp, 2.02_wp, 2.20_wp, &
                          3.48_wp, 2.83_wp, &
                          2.60_wp, &
                          2.37_wp, 2.43_wp, 2.40_wp, 2.21_wp, 2.43_wp, 2.44_wp, 2.45_wp, 2.44_wp, 2.45_wp, 2.45_wp, &
                          2.45_wp, 2.46_wp, 2.46_wp, 2.45_wp, &
                          [(2.45_wp, i=1, 9)], &
                          [(2.45_wp, i=1, 6)] &
                          ]

   !> D3 van-der-Waals radii
   real(wp), parameter :: d3_vdw_rad(1:94) = aatoau*[&
      & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &
      & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &
      & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &
      & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &
      & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &
      & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &
      & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &
      & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &
      & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &
      & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &
      & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &
      & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &
      & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &
      & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &
      & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &
      & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &
      & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &
      & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &
      & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &
      & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &
      & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &
      & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &
      & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &
      & 2.17695_wp, 2.21705_wp]

   !> Default value for unoptimized van-der-Waals radii
   real(wp), parameter :: cosmostub = 2.223_wp

   !> COSMO optimized van-der-Waals radii
   real(wp), parameter :: cosmo_vdw_rad(1:94) = aatoau*[ &
      & 1.3000_wp, 1.6380_wp, 1.5700_wp, 1.0530_wp, &   ! H-Be
      & 2.0480_wp, 2.0000_wp, 1.8300_wp, 1.7200_wp, &   ! B-O
      & 1.7200_wp, 1.8018_wp, 1.8000_wp, 1.6380_wp, &   ! F-Mg
      & 2.1530_wp, 2.2000_wp, 2.1060_wp, 2.1600_wp, &   ! Al-S
      & 2.0500_wp, 2.2000_wp, 2.2230_wp, cosmostub, &   ! Cl-Ca
      & cosmostub, 2.2930_wp, cosmostub, cosmostub, &   ! Sc-Cr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Mn-Ni
      & cosmostub, 1.6260_wp, cosmostub, 2.7000_wp, &   ! Cu-Ge
      & 2.3500_wp, 2.2000_wp, 2.1600_wp, 2.3630_wp, &   ! As-Kr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rb-Zr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Nb-Ru
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rh-Cd
      & 2.2580_wp, 2.5500_wp, 2.4100_wp, 2.4100_wp, &   ! In-Te
      & 2.3200_wp, 2.5270_wp, cosmostub, cosmostub, &   ! I-Ba
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! La-Nd
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Pm-Gd
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tb-Er
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tm-Hf
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ta-Os
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ir-Hg
      & cosmostub, 2.3600_wp, 2.4220_wp, 2.3050_wp, &   ! Tl-Po
      & 2.3630_wp, 2.5740_wp, cosmostub, cosmostub, &   ! At-Ra
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ac-U
      & cosmostub, cosmostub]                           ! Np-Pu

   !> Gaussian charge scheme radii (Bondi-based, uniformly scaled by 1.2)
   !> Ref: J. Phys. Chem. 2010, 133, 244111
   !> Base: Bondi radii (J. Phys. Chem. 1964, 68, 441-451) with H = 1.1 A,
   !> Mantina et al. (J. Phys. Chem. A 2009, 113, 5806-5812) for 16 missing
   !> main-group elements, and 2.0 A fallback for remaining elements.
   real(wp), parameter :: gauss_vdw_rad(max_elem_gauss) = (aatoau*1.2_wp)*[ &
      & 1.100_wp, 1.400_wp, &                                                       ! H -He
      & 1.820_wp, 1.530_wp, 1.920_wp, 1.700_wp, 1.550_wp, 1.520_wp, &              ! Li-O
      & 1.470_wp, 1.540_wp, &                                                       ! F -Ne
      & 2.270_wp, 1.730_wp, 1.840_wp, 2.100_wp, 1.800_wp, 1.800_wp, &              ! Na-S
      & 1.750_wp, 1.880_wp, &                                                       ! Cl-Ar
      & 2.750_wp, 2.310_wp, &                                                       ! K -Ca
      & [(2.000_wp, i=1, 7)], 1.630_wp, 1.400_wp, 1.390_wp, &                        ! Sc-Zn
      & 1.870_wp, 2.110_wp, 1.850_wp, 1.900_wp, 1.850_wp, 2.020_wp, &              ! Ga-Kr
      & 3.030_wp, 2.500_wp, &                                                       ! Rb-Sr
      & [(2.000_wp, i=1, 7)], 1.630_wp, 1.720_wp, 1.580_wp, &                        ! Y -Cd
      & 1.930_wp, 2.170_wp, 2.060_wp, 2.060_wp, 1.980_wp, 2.160_wp, &              ! In-Xe
      & 3.430_wp, 2.680_wp, &                                                       ! Cs-Ba
      & 2.000_wp, &                                                                 ! La
      & [(2.000_wp, i=1, 14)], &                                                     ! Ce-Lu
      & [(2.000_wp, i=1, 6)], 1.720_wp, 1.660_wp, 1.550_wp, &                        ! Hf-Hg
      & 1.960_wp, 2.020_wp, 2.070_wp, 1.970_wp, 2.020_wp, 2.200_wp, &              ! Tl-Rn
      & 3.480_wp, 2.830_wp, &                                                       ! Fr-Ra
      & 2.000_wp, &                                                                 ! Ac
      & 2.000_wp, 2.000_wp, 1.860_wp, [(2.000_wp, i=1, 11)], &                       ! Th-Lr
      & [(2.000_wp, i=1, 9)], &                                                      ! Rf-Cn
      & [(2.000_wp, i=1, 6)] &                                                       ! Nh-Og
   ]

   !> In case no van-der-Waals value is provided
   ! TODO: this is dangerous, better use some stub?
   real(wp), parameter :: missing = -1.0_wp

   real(wp), parameter :: bondi_vdw_rad(1:88) = aatoau*[ &
      & 1.10_wp, 1.40_wp, 1.81_wp, 1.53_wp, 1.92_wp, 1.70_wp, 1.55_wp, 1.52_wp, &  ! H-O
      & 1.47_wp, 1.54_wp, 2.27_wp, 1.73_wp, 1.84_wp, 2.10_wp, 1.80_wp, 1.80_wp, &  ! F-S
      & 1.75_wp, 1.88_wp, 2.75_wp, 2.31_wp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, missing, missing, missing, missing, missing, 1.87_wp, 2.11_wp, &  ! Mn-Ge
      & 1.85_wp, 1.90_wp, 1.83_wp, 2.02_wp, 3.03_wp, 2.49_wp, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 1.93_wp, 2.17_wp, 2.06_wp, 2.06_wp, 1.98_wp, 2.16_wp, 3.43_wp, 2.68_wp, &  ! I-Ba
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! La-Gd
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Tb-Hf
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Ta-Hg
      & 1.96_wp, 2.02_wp, 2.07_wp, 1.97_wp, 2.02_wp, 2.20_wp, 3.48_wp, 2.83_wp]    ! Tl-Ra

   !> Rahm, Hoffmann & Ashcroft (2016) atomic radii
   !> 0.001 e/bohr^3 isodensity surface of free atoms (DFT PBE0)
   !> Ref: Chem. Eur. J. 22, 14625 (2016)
   real(wp), parameter :: rahm_vdw_rad(1:96) = aatoau*[ &
      & 1.54_wp, 1.34_wp, 2.20_wp, 2.19_wp, 2.05_wp, 1.90_wp, 1.79_wp, 1.71_wp, &  ! H-O
      & 1.63_wp, 1.56_wp, 2.25_wp, 2.40_wp, 2.39_wp, 2.32_wp, 2.23_wp, 2.14_wp, &  ! F-S
      & 2.06_wp, 1.97_wp, 2.34_wp, 2.70_wp, 2.63_wp, 2.57_wp, 2.52_wp, 2.33_wp, &  ! Cl-Cr
      & 2.42_wp, 2.26_wp, 2.22_wp, 2.19_wp, 2.17_wp, 2.22_wp, 2.33_wp, 2.34_wp, &  ! Mn-Ge
      & 2.31_wp, 2.24_wp, 2.19_wp, 2.12_wp, 2.40_wp, 2.79_wp, 2.74_wp, 2.68_wp, &  ! As-Zr
      & 2.51_wp, 2.44_wp, 2.41_wp, 2.37_wp, 2.33_wp, 2.15_wp, 2.25_wp, 2.38_wp, &  ! Nb-Cd
      & 2.46_wp, 2.48_wp, 2.46_wp, 2.42_wp, 2.38_wp, 2.32_wp, 2.49_wp, 2.93_wp, &  ! In-Ba
      & 2.84_wp, 2.82_wp, 2.86_wp, 2.84_wp, 2.83_wp, 2.80_wp, 2.80_wp, 2.77_wp, &  ! La-Gd
      & 2.76_wp, 2.75_wp, 2.73_wp, 2.72_wp, 2.71_wp, 2.77_wp, 2.70_wp, 2.64_wp, &  ! Tb-Hf
      & 2.58_wp, 2.53_wp, 2.49_wp, 2.44_wp, 2.33_wp, 2.30_wp, 2.26_wp, 2.29_wp, &  ! Ta-Hg
      & 2.42_wp, 2.49_wp, 2.50_wp, 2.50_wp, 2.47_wp, 2.43_wp, 2.58_wp, 2.92_wp, &  ! Tl-Ra
      & 2.93_wp, 2.89_wp, 2.85_wp, 2.83_wp, 2.80_wp, 2.78_wp, 2.76_wp, 2.76_wp]    ! Ac-Cm

   public :: get_radius
   public :: get_radius_func

   interface get_radius
      module procedure get_radius_default
      module procedure get_atnum_int
      module procedure get_atnum_str
      module procedure get_atsym_int
      module procedure get_atsym_str
   end interface

   interface get_radius_func
      module procedure get_radius_func_default
      module procedure get_radius_func_int
      module procedure get_radius_func_str
   end interface

contains

   !> Look up the maximum valid atomic number for a given model tag.
   !> @param[in]  model   radius type integer tag
   !> @param[out] upper   upper bound on atomic number for that model
   pure subroutine get_upper_bound(model, upper)
      integer, intent(in)  :: model
      integer, intent(out) :: upper
      select case (model)
      case (rad_type%cpcm); upper = max_elem_cpcm
      case (rad_type%smd); upper = max_elem_smd
      case (rad_type%d3); upper = max_elem_d3
      case (rad_type%cosmo); upper = max_elem_cosmo
      case (rad_type%bondi); upper = max_elem_bondi
      case (rad_type%rahm); upper = max_elem_rahm
      case (rad_type%gauss); upper = max_elem_gauss
      case default; upper = -1
      end select
   end subroutine get_upper_bound

   !> Resolve a model name string to its integer tag.
   !> Returns -1 for unknown names.
   !> @param[in]  name   to_lower model name
   !> @param[out] tag    integer tag or -1
   subroutine resolve_model_name(name, tag)
      character(len=*), intent(in) :: name
      integer, intent(out) :: tag
      character(len=:), allocatable :: mstr

      mstr = adjustl(to_lower(trim(name)))
      select case (mstr)
      case ("cpcm"); tag = rad_type%cpcm
      case ("smd"); tag = rad_type%smd
      case ("d3"); tag = rad_type%d3
      case ("cosmo"); tag = rad_type%cosmo
      case ("bondi"); tag = rad_type%bondi
      case ("rahm"); tag = rad_type%rahm
      case ("gauss"); tag = rad_type%gauss
      case default; tag = -1
      end select
   end subroutine resolve_model_name

   !> Retrieve a radius from the raw data arrays for a validated
   !> (num, model) pair. No bounds checking is done here.
   !> @param[in]  num    atomic number (must be in range)
   !> @param[in]  model  radius type tag (must be valid)
   !> @param[out] rad    van-der-Waals radius
   pure subroutine fetch_radius(num, model, rad)
      integer, intent(in)   :: num
      integer, intent(in)   :: model
      real(wp), intent(out) :: rad

      select case (model)
      case (rad_type%smd); rad = smd_vdw_rad(num)
      case (rad_type%d3); rad = d3_vdw_rad(num)
      case (rad_type%cosmo); rad = cosmo_vdw_rad(num)
      case (rad_type%bondi); rad = bondi_vdw_rad(num)
      case (rad_type%rahm); rad = rahm_vdw_rad(num)
      case (rad_type%gauss); rad = gauss_vdw_rad(num)
      case default; rad = cpcm_vdw_rad(num)
      end select
   end subroutine fetch_radius

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius for an atomic number using the default
   !> (CPCM) radii set.
   !> @param[in]  num    atomic number
   !> @param[out] rad    radius in atomic units
   !> @param[out] error  error on invalid atomic number
   subroutine get_radius_default(num, rad, error)
      integer, intent(in) :: num
      real(wp), intent(out) :: rad
      type(error_type), allocatable, intent(out) :: error
      character(len=64) :: msg

      if (num < 1 .or. num > max_elem_cpcm) then
         rad = 0.0_wp
         write (msg, '(a,i0,a,i0,a)') &
            "Atomic number ", num, " out of range [1, ", max_elem_cpcm, "]"
         call fatal_error(error, trim(msg))
         return
      end if

      rad = cpcm_vdw_rad(num)

      if (rad < 0.0_wp) then
         write (msg, '(a,i0,a)') &
            "No valid CPCM radius for atomic number ", num, ""
         call fatal_error(error, trim(msg))
      end if
   end subroutine get_radius_default

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius for an atomic number and integer model tag.
   !> @param[in]  num    atomic number
   !> @param[in]  model  radius type tag (from rad_type)
   !> @param[out] rad    radius in atomic units
   !> @param[out] error  error on invalid number or model
   subroutine get_atnum_int(num, model, rad, error)
      integer, intent(in)  :: num
      integer, intent(in)  :: model
      real(wp), intent(out) :: rad
      type(error_type), allocatable, intent(out) :: error
      integer :: upper
      character(len=128) :: msg

      call get_upper_bound(model, upper)
      if (upper < 0) then
         rad = 0.0_wp
         write (msg, '(a,i0)') "Unknown radius type: ", model
         call fatal_error(error, trim(msg))
         return
      end if

      if (num < 1 .or. num > upper) then
         rad = 0.0_wp
         write (msg, '(a,i0,a,i0,a)') &
            "Atomic number ", num, " out of range [1, ", upper, "]"
         call fatal_error(error, trim(msg))
         return
      end if

      call fetch_radius(num, model, rad)

      if (rad < 0.0_wp) then
         write (msg, '(a,i0)') &
            "No valid radius for atomic number ", num
         call fatal_error(error, trim(msg))
      end if
   end subroutine get_atnum_int

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius for an atomic number and model name string.
   !> @param[in]  num        atomic number
   !> @param[in]  model_name radius set name (cpcm, smd, d3, cosmo, bondi)
   !> @param[out] rad        radius in atomic units
   !> @param[out] error      error on invalid number or name
   subroutine get_atnum_str(num, model_name, rad, error)
      integer, intent(in)          :: num
      character(len=*), intent(in) :: model_name
      real(wp), intent(out) :: rad
      type(error_type), allocatable, intent(out) :: error
      integer :: tag
      character(len=128) :: msg

      call resolve_model_name(model_name, tag)
      if (tag < 0) then
         rad = 0.0_wp
         write (msg, '(a,a,a)') &
            "Unknown radius type: '", trim(model_name), "'"
         call fatal_error(error, trim(msg))
         return
      end if

      call get_atnum_int(num, tag, rad, error)
   end subroutine get_atnum_str

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius for an element symbol and integer model tag.
   !> @param[in]  sym    element symbol (e.g. 'H', 'He')
   !> @param[in]  model  radius type tag (from rad_type)
   !> @param[out] rad    radius in atomic units
   !> @param[out] error  error on invalid symbol or model
   subroutine get_atsym_int(sym, model, rad, error)
      character(len=*), intent(in) :: sym
      integer, intent(in)          :: model
      real(wp), intent(out)        :: rad
      type(error_type), allocatable, intent(out) :: error
      integer :: num
      character(len=128) :: msg

      call symbol_to_number(num, trim(sym))
      if (num < 1) then
         rad = 0.0_wp
         write (msg, '(a,a,a)') &
            "Unknown element symbol: '", trim(sym), "'"
         call fatal_error(error, trim(msg))
         return
      end if

      call get_atnum_int(num, model, rad, error)
   end subroutine get_atsym_int

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius for an element symbol and model name string.
   !> @param[in]  sym        element symbol (e.g. 'H', 'He')
   !> @param[in]  model_name radius set name (cpcm, smd, d3, cosmo, bondi)
   !> @param[out] rad        radius in atomic units
   !> @param[out] error      error on invalid symbol or name
   subroutine get_atsym_str(sym, model_name, rad, error)
      character(len=*), intent(in) :: sym
      character(len=*), intent(in) :: model_name
      real(wp), intent(out)        :: rad
      type(error_type), allocatable, intent(out) :: error
      integer :: num, tag
      character(len=128) :: msg

      call resolve_model_name(model_name, tag)
      if (tag < 0) then
         rad = 0.0_wp
         write (msg, '(a,a,a)') &
            "Unknown radius type: '", trim(model_name), "'"
         call fatal_error(error, trim(msg))
         return
      end if

      call symbol_to_number(num, trim(sym))
      if (num < 1) then
         rad = 0.0_wp
         write (msg, '(a,a,a)') &
            "Unknown element symbol: '", trim(sym), "'"
         call fatal_error(error, trim(msg))
         return
      end if

      call get_atnum_int(num, tag, rad, error)
   end subroutine get_atsym_str

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius as a function (for test/debug use).
   !> Uses the default (CPCM) radii set. Stops on error.
   !> @param[in]  num    atomic number
   !> @return     radius in atomic units
   function get_radius_func_default(num) result(rad)
      !> Atomic number
      integer, intent(in) :: num
      !> Radius in atomic units
      real(wp) :: rad
      !> Error handler (local)
      type(error_type), allocatable :: error

      call get_radius_default(num, rad, error)
      if (allocated(error)) then
         rad = missing
      end if
   end function get_radius_func_default

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius as a function with model tag (for test/debug use).
   !> Stops on error.
   !> @param[in]  num    atomic number
   !> @param[in]  model  radius type tag (from rad_type)
   !> @return     radius in atomic units
   function get_radius_func_int(num, model) result(rad)
      !> Atomic number
      integer, intent(in) :: num
      !> Radius model tag
      integer, intent(in) :: model
      !> Radius in atomic units
      real(wp) :: rad
      !> Error handler (local)
      type(error_type), allocatable :: error

      call get_atnum_int(num, model, rad, error)
      if (allocated(error)) then
         rad = missing
      end if
   end function get_radius_func_int

   !--------------------------------------------------------------------
   !> Get van-der-Waals radius as a function with model name (for test/debug use).
   !> Stops on error.
   !> @param[in]  num         atomic number
   !> @param[in]  model_name  radius set name (cpcm, smd, d3, cosmo, bondi)
   !> @return     radius in atomic units
   function get_radius_func_str(num, model_name) result(rad)
      !> Atomic number
      integer, intent(in) :: num
      !> Radius model name
      character(len=*), intent(in) :: model_name
      !> Radius in atomic units
      real(wp) :: rad
      !> Error handler (local)
      type(error_type), allocatable :: error

      call get_atnum_str(num, model_name, rad, error)
      if (allocated(error)) then
         rad = missing
      end if
   end function get_radius_func_str

end module moist_data_radii_legacy
