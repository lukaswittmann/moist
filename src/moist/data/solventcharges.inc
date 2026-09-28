! Solvent atomic charges and MBIS Cartesian multipoles
! Source: solvent_data_final.zip, dataset.json schema_version 1.0
!
! All single-point calculations:
! wB97M-V, def2-TZVPP/def2/J, RIJCOSX, TightSCF, DEFGRID3
! IDs 22, 23, 24, 25, 106, 129, 130, 131: SF-X2C-1e/x2c-TZVPPall, x2c/J
!
! Environments:
! - gas (no CPCM),
! - solvent (CPCM with solvent's dielectric),
! - conductor (CPCM infinite dielectric)
!
! Charges (in e):
! - RESP is ORCA hyperbolically restrained CHELPG (a=0.0005, b=0.1);
! - CHELPG uses the same 0.3 Angstrom grid and 2.8 Angstrom rmax
!
subroutine get_solvent_charges(solvent_id, environment, model, charges, error)
   integer, intent(in) :: solvent_id
   character(len=*), intent(in) :: environment
   character(len=*), intent(in) :: model
   real(wp), allocatable, intent(out) :: charges(:)
   type(error_type), allocatable, intent(out) :: error
   character(len=64) :: msg

   select case (solvent_id)
   case (1)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073789_wp, 0.15738_wp, -0.073736_wp, -0.073833_wp, -0.071001_wp, &
                      0.044998_wp, 0.044985_wp, 0.044996_wp]
         case ("resp")
            charges = [ &
                      -0.015318_wp, -0.065822_wp, -0.015794_wp, -0.015585_wp, -0.478797_wp, &
                      0.197087_wp, 0.196492_wp, 0.197736_wp]
         case ("mbis")
            charges = [ &
                      -0.110492_wp, 0.237688_wp, -0.11041_wp, -0.110549_wp, -0.478119_wp, &
                      0.190635_wp, 0.190618_wp, 0.190628_wp]
         case ("chelpg")
            charges = [ &
                      -0.015089_wp, -0.063728_wp, -0.01556_wp, -0.015359_wp, -0.493899_wp, &
                      0.201195_wp, 0.200601_wp, 0.201839_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.082124_wp, 0.159064_wp, -0.082081_wp, -0.082187_wp, -0.068378_wp, &
                      0.051881_wp, 0.051907_wp, 0.051917_wp]
         case ("resp")
            charges = [ &
                      -0.031026_wp, -0.072655_wp, -0.031486_wp, -0.031332_wp, -0.409442_wp, &
                      0.191903_wp, 0.19137_wp, 0.192667_wp]
         case ("mbis")
            charges = [ &
                      -0.122228_wp, 0.249567_wp, -0.122162_wp, -0.122312_wp, -0.490998_wp, &
                      0.202677_wp, 0.202723_wp, 0.202731_wp]
         case ("chelpg")
            charges = [ &
                      -0.030828_wp, -0.070465_wp, -0.031283_wp, -0.031137_wp, -0.424395_wp, &
                      0.19596_wp, 0.19543_wp, 0.196719_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.083861_wp, 0.159417_wp, -0.083821_wp, -0.083929_wp, -0.067709_wp, &
                      0.053274_wp, 0.053309_wp, 0.053319_wp]
         case ("resp")
            charges = [ &
                      -0.034551_wp, -0.074242_wp, -0.035008_wp, -0.034866_wp, -0.391377_wp, &
                      0.189924_wp, 0.189405_wp, 0.190714_wp]
         case ("mbis")
            charges = [ &
                      -0.124646_wp, 0.251827_wp, -0.124585_wp, -0.124738_wp, -0.492945_wp, &
                      0.204985_wp, 0.205047_wp, 0.205054_wp]
         case ("chelpg")
            charges = [ &
                      -0.034361_wp, -0.072036_wp, -0.034813_wp, -0.034679_wp, -0.406279_wp, &
                      0.193966_wp, 0.19345_wp, 0.194752_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (2)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.081782_wp, 0.08894_wp, -0.080507_wp, 0.014673_wp, -0.10339_wp, &
                      0.056975_wp, 0.051936_wp, 0.053156_wp]
         case ("resp")
            charges = [ &
                      -0.073256_wp, -0.144063_wp, -0.05175_wp, -0.269254_wp, -0.103881_wp, &
                      0.233919_wp, 0.192107_wp, 0.216177_wp]
         case ("mbis")
            charges = [ &
                      -0.147507_wp, 0.04846_wp, -0.126938_wp, -0.181779_wp, -0.16476_wp, &
                      0.19443_wp, 0.185745_wp, 0.192349_wp]
         case ("chelpg")
            charges = [ &
                      -0.072923_wp, -0.145056_wp, -0.051219_wp, -0.275634_wp, -0.102776_wp, &
                      0.234961_wp, 0.194207_wp, 0.21844_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.08906_wp, 0.09449_wp, -0.098733_wp, 0.020836_wp, -0.122957_wp, &
                      0.06972_wp, 0.064376_wp, 0.061327_wp]
         case ("resp")
            charges = [ &
                      -0.092018_wp, -0.130504_wp, -0.087458_wp, -0.228378_wp, -0.144782_wp, &
                      0.256523_wp, 0.20854_wp, 0.218078_wp]
         case ("mbis")
            charges = [ &
                      -0.159464_wp, 0.059279_wp, -0.15268_wp, -0.178933_wp, -0.191571_wp, &
                      0.213929_wp, 0.204852_wp, 0.204589_wp]
         case ("chelpg")
            charges = [ &
                      -0.091733_wp, -0.131287_wp, -0.086985_wp, -0.234656_wp, -0.143704_wp, &
                      0.257482_wp, 0.210591_wp, 0.22029_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.090198_wp, 0.095657_wp, -0.102702_wp, 0.022178_wp, -0.126917_wp, &
                      0.072181_wp, 0.066998_wp, 0.062803_wp]
         case ("resp")
            charges = [ &
                      -0.095424_wp, -0.127624_wp, -0.095396_wp, -0.218781_wp, -0.153362_wp, &
                      0.260924_wp, 0.211954_wp, 0.217709_wp]
         case ("mbis")
            charges = [ &
                      -0.161454_wp, 0.061597_wp, -0.15826_wp, -0.177951_wp, -0.197028_wp, &
                      0.217619_wp, 0.208798_wp, 0.20668_wp]
         case ("chelpg")
            charges = [ &
                      -0.095146_wp, -0.128374_wp, -0.094934_wp, -0.225025_wp, -0.152293_wp, &
                      0.26187_wp, 0.213995_wp, 0.219908_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (3)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.050462_wp, 0.002338_wp, -0.002651_wp, -0.044898_wp, -0.049881_wp, &
                      0.002735_wp, -0.068388_wp, -0.068546_wp, -0.068929_wp, 0.028412_wp, &
                      0.026823_wp, 0.029657_wp, 0.029869_wp, 0.029157_wp, 0.029119_wp, &
                      0.029487_wp, 0.029533_wp, 0.028453_wp, 0.029149_wp, 0.02945_wp, &
                      0.02958_wp]
         case ("resp")
            charges = [ &
                      -0.509833_wp, 0.260459_wp, 0.163662_wp, -0.258934_wp, -0.358161_wp, &
                      0.396521_wp, -0.461564_wp, -0.44642_wp, -0.453138_wp, 0.125219_wp, &
                      0.203232_wp, 0.172906_wp, 0.179215_wp, 0.123897_wp, 0.122038_wp, &
                      0.123875_wp, 0.124798_wp, 0.126268_wp, 0.122059_wp, 0.118818_wp, &
                      0.125083_wp]
         case ("mbis")
            charges = [ &
                      -0.255874_wp, 0.090841_wp, 0.060788_wp, -0.176348_wp, -0.21206_wp, &
                      0.147853_wp, -0.456806_wp, -0.435496_wp, -0.430336_wp, 0.137175_wp, &
                      0.136725_wp, 0.132006_wp, 0.13767_wp, 0.138635_wp, 0.138872_wp, &
                      0.140443_wp, 0.139916_wp, 0.136889_wp, 0.142641_wp, 0.143112_wp, &
                      0.14336_wp]
         case ("chelpg")
            charges = [ &
                      -0.533535_wp, 0.2806_wp, 0.158354_wp, -0.248641_wp, -0.377323_wp, &
                      0.42402_wp, -0.484642_wp, -0.463423_wp, -0.468663_wp, 0.129275_wp, &
                      0.207947_wp, 0.171945_wp, 0.183224_wp, 0.127381_wp, 0.125513_wp, &
                      0.128046_wp, 0.128994_wp, 0.130299_wp, 0.127104_wp, 0.12356_wp, &
                      0.129966_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.055135_wp, -0.000713_wp, -0.006017_wp, -0.050399_wp, -0.055706_wp, &
                      -0.001404_wp, -0.069641_wp, -0.069644_wp, -0.070123_wp, 0.030711_wp, &
                      0.030754_wp, 0.033175_wp, 0.033114_wp, 0.031542_wp, 0.03148_wp, &
                      0.031646_wp, 0.031574_wp, 0.030754_wp, 0.031287_wp, 0.031548_wp, &
                      0.031203_wp]
         case ("resp")
            charges = [ &
                      -0.526951_wp, 0.259688_wp, 0.166545_wp, -0.277285_wp, -0.37501_wp, &
                      0.395696_wp, -0.47768_wp, -0.465843_wp, -0.474368_wp, 0.133999_wp, &
                      0.216698_wp, 0.185387_wp, 0.191123_wp, 0.13262_wp, 0.130588_wp, &
                      0.130947_wp, 0.132291_wp, 0.13492_wp, 0.128756_wp, 0.12574_wp, &
                      0.13214_wp]
         case ("mbis")
            charges = [ &
                      -0.272875_wp, 0.090925_wp, 0.05969_wp, -0.191101_wp, -0.22787_wp, &
                      0.148651_wp, -0.467839_wp, -0.447876_wp, -0.442612_wp, 0.1432_wp, &
                      0.147988_wp, 0.142566_wp, 0.147872_wp, 0.144772_wp, 0.144966_wp, &
                      0.146208_wp, 0.14553_wp, 0.142916_wp, 0.148239_wp, 0.148677_wp, &
                      0.147979_wp]
         case ("chelpg")
            charges = [ &
                      -0.55064_wp, 0.279825_wp, 0.161345_wp, -0.267057_wp, -0.394159_wp, &
                      0.423223_wp, -0.500891_wp, -0.482967_wp, -0.490063_wp, 0.138094_wp, &
                      0.221411_wp, 0.184442_wp, 0.195132_wp, 0.136136_wp, 0.134094_wp, &
                      0.135151_wp, 0.136527_wp, 0.13899_wp, 0.133834_wp, 0.130515_wp, &
                      0.137057_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.059108_wp, -0.003346_wp, -0.008996_wp, -0.055251_wp, -0.060947_wp, &
                      -0.005187_wp, -0.07062_wp, -0.070388_wp, -0.070988_wp, 0.03283_wp, &
                      0.034171_wp, 0.036042_wp, 0.035681_wp, 0.033725_wp, 0.033636_wp, &
                      0.033539_wp, 0.033272_wp, 0.032863_wp, 0.033198_wp, 0.033301_wp, &
                      0.032582_wp]
         case ("resp")
            charges = [ &
                      -0.538652_wp, 0.256719_wp, 0.167468_wp, -0.29232_wp, -0.388455_wp, &
                      0.391521_wp, -0.489443_wp, -0.480732_wp, -0.490671_wp, 0.141559_wp, &
                      0.227647_wp, 0.195555_wp, 0.200736_wp, 0.14017_wp, 0.137937_wp, &
                      0.136738_wp, 0.13827_wp, 0.14231_wp, 0.13434_wp, 0.131312_wp, &
                      0.137991_wp]
         case ("mbis")
            charges = [ &
                      -0.286989_wp, 0.090658_wp, 0.058309_wp, -0.203647_wp, -0.241624_wp, &
                      0.148486_wp, -0.476454_wp, -0.457768_wp, -0.452446_wp, 0.148479_wp, &
                      0.15761_wp, 0.151355_wp, 0.156315_wp, 0.150106_wp, 0.150251_wp, &
                      0.15105_wp, 0.150099_wp, 0.148177_wp, 0.153037_wp, 0.15323_wp, &
                      0.151775_wp]
         case ("chelpg")
            charges = [ &
                      -0.562319_wp, 0.276846_wp, 0.162347_wp, -0.28214_wp, -0.407585_wp, &
                      0.41906_wp, -0.51275_wp, -0.497952_wp, -0.506491_wp, 0.145683_wp, &
                      0.232358_wp, 0.19462_wp, 0.204745_wp, 0.14371_wp, 0.141466_wp, &
                      0.140969_wp, 0.142536_wp, 0.14641_wp, 0.139443_wp, 0.136112_wp, &
                      0.142932_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (4)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.011151_wp, 0.049173_wp, -0.10947_wp, 0.049175_wp, 0.011135_wp, &
                      0.049169_wp, -0.109482_wp, 0.049163_wp]
         case ("resp")
            charges = [ &
                      -0.166458_wp, 0.157944_wp, -0.147542_wp, 0.157539_wp, -0.175608_wp, &
                      0.159975_wp, -0.145255_wp, 0.159405_wp]
         case ("mbis")
            charges = [ &
                      -0.201559_wp, 0.184093_wp, -0.166578_wp, 0.184074_wp, -0.201604_wp, &
                      0.184079_wp, -0.166581_wp, 0.184088_wp]
         case ("chelpg")
            charges = [ &
                      -0.170236_wp, 0.159548_wp, -0.146967_wp, 0.159137_wp, -0.179361_wp, &
                      0.161567_wp, -0.144681_wp, 0.160995_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.017138_wp, 0.056272_wp, -0.129658_wp, 0.056278_wp, 0.017121_wp, &
                      0.056271_wp, -0.129669_wp, 0.056261_wp]
         case ("resp")
            charges = [ &
                      -0.130967_wp, 0.160843_wp, -0.189048_wp, 0.160524_wp, -0.139461_wp, &
                      0.162768_wp, -0.186849_wp, 0.16219_wp]
         case ("mbis")
            charges = [ &
                      -0.193588_wp, 0.194038_wp, -0.194442_wp, 0.194024_wp, -0.19363_wp, &
                      0.194027_wp, -0.194445_wp, 0.19403_wp]
         case ("chelpg")
            charges = [ &
                      -0.134535_wp, 0.162367_wp, -0.188514_wp, 0.162042_wp, -0.143054_wp, &
                      0.164292_wp, -0.18631_wp, 0.163712_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.019182_wp, 0.058449_wp, -0.136056_wp, 0.058456_wp, 0.019164_wp, &
                      0.058448_wp, -0.136066_wp, 0.058437_wp]
         case ("resp")
            charges = [ &
                      -0.117606_wp, 0.161084_wp, -0.202971_wp, 0.160797_wp, -0.125854_wp, &
                      0.16297_wp, -0.200807_wp, 0.162389_wp]
         case ("mbis")
            charges = [ &
                      -0.190517_wp, 0.19694_wp, -0.20332_wp, 0.196928_wp, -0.190558_wp, &
                      0.196931_wp, -0.203323_wp, 0.196931_wp]
         case ("chelpg")
            charges = [ &
                      -0.121082_wp, 0.162572_wp, -0.202452_wp, 0.162279_wp, -0.129382_wp, &
                      0.164464_wp, -0.200279_wp, 0.163882_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (5)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.021879_wp, 0.021879_wp, -0.12039_wp, -0.120388_wp, 0.049478_wp, &
                      0.049033_wp, 0.049033_wp, 0.049477_wp]
         case ("resp")
            charges = [ &
                      -0.119421_wp, -0.140299_wp, -0.152679_wp, -0.14908_wp, 0.128856_wp, &
                      0.146825_wp, 0.152195_wp, 0.133602_wp]
         case ("mbis")
            charges = [ &
                      -0.145495_wp, -0.14549_wp, -0.185784_wp, -0.185781_wp, 0.159352_wp, &
                      0.171928_wp, 0.171923_wp, 0.159348_wp]
         case ("chelpg")
            charges = [ &
                      -0.12246_wp, -0.143616_wp, -0.152041_wp, -0.148403_wp, 0.12999_wp, &
                      0.148145_wp, 0.153579_wp, 0.134805_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.030374_wp, 0.030374_wp, -0.153165_wp, -0.153163_wp, 0.065217_wp, &
                      0.057573_wp, 0.057574_wp, 0.065216_wp]
         case ("resp")
            charges = [ &
                      -0.097947_wp, -0.118116_wp, -0.209994_wp, -0.206423_wp, 0.154479_wp, &
                      0.156886_wp, 0.162093_wp, 0.159022_wp]
         case ("mbis")
            charges = [ &
                      -0.137277_wp, -0.137274_wp, -0.229927_wp, -0.229924_wp, 0.18331_wp, &
                      0.183895_wp, 0.183892_wp, 0.183306_wp]
         case ("chelpg")
            charges = [ &
                      -0.100809_wp, -0.121372_wp, -0.209392_wp, -0.205765_wp, 0.155556_wp, &
                      0.158143_wp, 0.163442_wp, 0.160196_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.031677_wp, 0.031676_wp, -0.157909_wp, -0.157906_wp, 0.067538_wp, &
                      0.058693_wp, 0.058694_wp, 0.067538_wp]
         case ("resp")
            charges = [ &
                      -0.094683_wp, -0.114734_wp, -0.218437_wp, -0.21487_wp, 0.158385_wp, &
                      0.158133_wp, 0.163312_wp, 0.162894_wp]
         case ("mbis")
            charges = [ &
                      -0.135816_wp, -0.135814_wp, -0.236342_wp, -0.236338_wp, 0.186786_wp, &
                      0.185374_wp, 0.18537_wp, 0.186782_wp]
         case ("chelpg")
            charges = [ &
                      -0.097505_wp, -0.117972_wp, -0.217844_wp, -0.214218_wp, 0.159449_wp, &
                      0.159377_wp, 0.164653_wp, 0.16406_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (6)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.040851_wp, 0.037771_wp, -0.182889_wp, 0.037727_wp, 0.040888_wp, &
                      -0.250254_wp, -0.250066_wp, 0.039437_wp, 0.025019_wp, 0.026656_wp, &
                      0.027943_wp, 0.027903_wp, 0.026666_wp, 0.025075_wp, 0.039465_wp, &
                      0.143724_wp, 0.144083_wp]
         case ("resp")
            charges = [ &
                      0.109479_wp, -0.114958_wp, -0.200981_wp, -0.110131_wp, 0.106772_wp, &
                      -0.561356_wp, -0.560358_wp, 0.092866_wp, 0.038561_wp, 0.105572_wp, &
                      0.070792_wp, 0.070023_wp, 0.104777_wp, 0.039214_wp, 0.092718_wp, &
                      0.358637_wp, 0.358373_wp]
         case ("mbis")
            charges = [ &
                      0.028854_wp, -0.039399_wp, -0.365497_wp, -0.039051_wp, 0.028596_wp, &
                      -0.608034_wp, -0.607901_wp, 0.108568_wp, 0.078754_wp, 0.094688_wp, &
                      0.08876_wp, 0.088585_wp, 0.094767_wp, 0.078791_wp, 0.108683_wp, &
                      0.43044_wp, 0.430395_wp]
         case ("chelpg")
            charges = [ &
                      0.115836_wp, -0.122993_wp, -0.198547_wp, -0.11797_wp, 0.112914_wp, &
                      -0.563087_wp, -0.562042_wp, 0.091681_wp, 0.037151_wp, 0.107686_wp, &
                      0.072833_wp, 0.072017_wp, 0.106857_wp, 0.037856_wp, 0.091581_wp, &
                      0.35925_wp, 0.358976_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.03894_wp, 0.039885_wp, -0.188798_wp, 0.039839_wp, 0.039005_wp, &
                      -0.277965_wp, -0.277763_wp, 0.043845_wp, 0.029119_wp, 0.030872_wp, &
                      0.037475_wp, 0.03747_wp, 0.030792_wp, 0.029252_wp, 0.043847_wp, &
                      0.15188_wp, 0.152304_wp]
         case ("resp")
            charges = [ &
                      0.118321_wp, -0.109036_wp, -0.227889_wp, -0.104462_wp, 0.116024_wp, &
                      -0.642824_wp, -0.642052_wp, 0.100664_wp, 0.047136_wp, 0.113022_wp, &
                      0.087352_wp, 0.0867_wp, 0.112067_wp, 0.047889_wp, 0.100389_wp, &
                      0.398404_wp, 0.398295_wp]
         case ("mbis")
            charges = [ &
                      0.028234_wp, -0.042886_wp, -0.378912_wp, -0.042507_wp, 0.028039_wp, &
                      -0.67038_wp, -0.670329_wp, 0.117815_wp, 0.087065_wp, 0.101853_wp, &
                      0.105589_wp, 0.105479_wp, 0.101771_wp, 0.087216_wp, 0.117866_wp, &
                      0.462005_wp, 0.462081_wp]
         case ("chelpg")
            charges = [ &
                      0.12469_wp, -0.117133_wp, -0.225438_wp, -0.11237_wp, 0.122194_wp, &
                      -0.64456_wp, -0.643745_wp, 0.09948_wp, 0.045728_wp, 0.115154_wp, &
                      0.089413_wp, 0.088716_wp, 0.114165_wp, 0.046529_wp, 0.099248_wp, &
                      0.399024_wp, 0.398905_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.038886_wp, 0.039964_wp, -0.188971_wp, 0.039918_wp, 0.038951_wp, &
                      -0.278792_wp, -0.27859_wp, 0.04396_wp, 0.029236_wp, 0.031003_wp, &
                      0.037785_wp, 0.037782_wp, 0.03092_wp, 0.029371_wp, 0.043961_wp, &
                      0.152094_wp, 0.15252_wp]
         case ("resp")
            charges = [ &
                      0.118551_wp, -0.108706_wp, -0.228753_wp, -0.104142_wp, 0.116266_wp, &
                      -0.645267_wp, -0.644501_wp, 0.100871_wp, 0.047392_wp, 0.113217_wp, &
                      0.087864_wp, 0.087217_wp, 0.112256_wp, 0.04815_wp, 0.100592_wp, &
                      0.39955_wp, 0.399445_wp]
         case ("mbis")
            charges = [ &
                      0.028236_wp, -0.042952_wp, -0.379319_wp, -0.042572_wp, 0.028043_wp, &
                      -0.672223_wp, -0.672175_wp, 0.118059_wp, 0.087301_wp, 0.102064_wp, &
                      0.106125_wp, 0.106017_wp, 0.101976_wp, 0.087456_wp, 0.118108_wp, &
                      0.462888_wp, 0.462966_wp]
         case ("chelpg")
            charges = [ &
                      0.124913_wp, -0.116797_wp, -0.226305_wp, -0.112045_wp, 0.122429_wp, &
                      -0.647004_wp, -0.646196_wp, 0.099688_wp, 0.045986_wp, 0.115348_wp, &
                      0.089925_wp, 0.089234_wp, 0.114353_wp, 0.046791_wp, 0.099452_wp, &
                      0.400171_wp, 0.400057_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (7)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.199328_wp, 0.036076_wp, 0.036074_wp, 0.024499_wp, 0.024494_wp, &
                      0.039088_wp, 0.039081_wp, -0.199325_wp, 0.036081_wp, 0.036079_wp, &
                      0.024501_wp, 0.024499_wp, 0.039081_wp, 0.039089_wp]
         case ("resp")
            charges = [ &
                      -0.357966_wp, 0.033923_wp, 0.045804_wp, 0.073522_wp, 0.070565_wp, &
                      0.068819_wp, 0.06699_wp, -0.364542_wp, 0.052562_wp, 0.039994_wp, &
                      0.068503_wp, 0.071171_wp, 0.064537_wp, 0.066119_wp]
         case ("mbis")
            charges = [ &
                      -0.354971_wp, -0.013414_wp, -0.013353_wp, 0.0815_wp, 0.08149_wp, &
                      0.109383_wp, 0.109364_wp, -0.354989_wp, -0.01337_wp, -0.013399_wp, &
                      0.081501_wp, 0.081493_wp, 0.109369_wp, 0.109384_wp]
         case ("chelpg")
            charges = [ &
                      -0.358191_wp, 0.033781_wp, 0.046952_wp, 0.073516_wp, 0.070258_wp, &
                      0.068803_wp, 0.06673_wp, -0.36533_wp, 0.054041_wp, 0.040167_wp, &
                      0.068083_wp, 0.071057_wp, 0.064152_wp, 0.065981_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.210842_wp, 0.035742_wp, 0.035738_wp, 0.027207_wp, 0.027196_wp, &
                      0.042474_wp, 0.042466_wp, -0.210836_wp, 0.035746_wp, 0.035746_wp, &
                      0.027204_wp, 0.027207_wp, 0.042466_wp, 0.042475_wp]
         case ("resp")
            charges = [ &
                      -0.394434_wp, 0.043984_wp, 0.056312_wp, 0.077194_wp, 0.074115_wp, &
                      0.073234_wp, 0.071323_wp, -0.401211_wp, 0.06323_wp, 0.050096_wp, &
                      0.072016_wp, 0.074837_wp, 0.068801_wp, 0.070502_wp]
         case ("mbis")
            charges = [ &
                      -0.375455_wp, -0.015102_wp, -0.01504_wp, 0.086836_wp, 0.086814_wp, &
                      0.115981_wp, 0.115962_wp, -0.37547_wp, -0.015055_wp, -0.015086_wp, &
                      0.086826_wp, 0.086829_wp, 0.115966_wp, 0.115982_wp]
         case ("chelpg")
            charges = [ &
                      -0.394771_wp, 0.04407_wp, 0.0576_wp, 0.077107_wp, 0.073747_wp, &
                      0.073146_wp, 0.071008_wp, -0.402072_wp, 0.064829_wp, 0.050476_wp, &
                      0.071542_wp, 0.074649_wp, 0.068368_wp, 0.070301_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.221846_wp, 0.035485_wp, 0.035479_wp, 0.029747_wp, 0.029728_wp, &
                      0.045696_wp, 0.045688_wp, -0.221837_wp, 0.035488_wp, 0.035489_wp, &
                      0.029737_wp, 0.029748_wp, 0.045688_wp, 0.045698_wp]
         case ("resp")
            charges = [ &
                      -0.430455_wp, 0.054648_wp, 0.067393_wp, 0.080342_wp, 0.077147_wp, &
                      0.077354_wp, 0.075367_wp, -0.437422_wp, 0.07448_wp, 0.060774_wp, &
                      0.075008_wp, 0.077991_wp, 0.072774_wp, 0.074598_wp]
         case ("mbis")
            charges = [ &
                      -0.395429_wp, -0.016233_wp, -0.01617_wp, 0.091758_wp, 0.091723_wp, &
                      0.12218_wp, 0.122163_wp, -0.395438_wp, -0.016185_wp, -0.016218_wp, &
                      0.091737_wp, 0.091752_wp, 0.122166_wp, 0.122182_wp]
         case ("chelpg")
            charges = [ &
                      -0.430901_wp, 0.054958_wp, 0.068799_wp, 0.080178_wp, 0.076727_wp, &
                      0.077196_wp, 0.075002_wp, -0.438348_wp, 0.076176_wp, 0.061352_wp, &
                      0.07449_wp, 0.077735_wp, 0.072301_wp, 0.074337_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (8)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.067243_wp, -0.002496_wp, -0.072745_wp, 0.011646_wp, -0.130271_wp, &
                      0.025274_wp, 0.026439_wp, 0.02959_wp, 0.027578_wp, 0.02399_wp, &
                      0.022651_wp, 0.026807_wp, 0.038927_wp, 0.039863_wp]
         case ("resp")
            charges = [ &
                      -0.531998_wp, 0.571382_wp, -0.516143_wp, -0.451817_wp, -0.162705_wp, &
                      0.125816_wp, 0.119501_wp, 0.127851_wp, -0.028351_wp, 0.13197_wp, &
                      0.120898_wp, 0.117275_wp, 0.192569_wp, 0.183752_wp]
         case ("mbis")
            charges = [ &
                      -0.459437_wp, 0.030446_wp, -0.447703_wp, -0.216957_wp, -0.205324_wp, &
                      0.139902_wp, 0.137947_wp, 0.14603_wp, 0.120415_wp, 0.134089_wp, &
                      0.159564_wp, 0.136532_wp, 0.161886_wp, 0.162617_wp]
         case ("chelpg")
            charges = [ &
                      -0.551952_wp, 0.588063_wp, -0.535005_wp, -0.46425_wp, -0.161385_wp, &
                      0.130801_wp, 0.124055_wp, 0.132636_wp, -0.029729_wp, 0.136578_wp, &
                      0.125468_wp, 0.121617_wp, 0.196082_wp, 0.18702_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.066207_wp, -0.000646_wp, -0.074762_wp, 0.018452_wp, -0.16524_wp, &
                      0.028474_wp, 0.029847_wp, 0.031657_wp, 0.029066_wp, 0.026795_wp, &
                      0.017506_wp, 0.029377_wp, 0.047399_wp, 0.048291_wp]
         case ("resp")
            charges = [ &
                      -0.571882_wp, 0.619448_wp, -0.548097_wp, -0.415856_wp, -0.224659_wp, &
                      0.137092_wp, 0.129529_wp, 0.137878_wp, -0.036713_wp, 0.140646_wp, &
                      0.117431_wp, 0.12653_wp, 0.199415_wp, 0.18924_wp]
         case ("mbis")
            charges = [ &
                      -0.469844_wp, 0.036324_wp, -0.46114_wp, -0.205311_wp, -0.251796_wp, &
                      0.146254_wp, 0.144558_wp, 0.150562_wp, 0.122188_wp, 0.141528_wp, &
                      0.153231_wp, 0.143751_wp, 0.174594_wp, 0.175108_wp]
         case ("chelpg")
            charges = [ &
                      -0.592063_wp, 0.636294_wp, -0.567073_wp, -0.428242_wp, -0.223353_wp, &
                      0.142127_wp, 0.134135_wp, 0.142717_wp, -0.038132_wp, 0.145275_wp, &
                      0.122021_wp, 0.130899_wp, 0.202907_wp, 0.192487_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.065879_wp, -0.000284_wp, -0.075108_wp, 0.019885_wp, -0.172083_wp, &
                      0.029173_wp, 0.030543_wp, 0.032007_wp, 0.029245_wp, 0.02738_wp, &
                      0.016321_wp, 0.029843_wp, 0.049042_wp, 0.049924_wp]
         case ("resp")
            charges = [ &
                      -0.578255_wp, 0.626122_wp, -0.552732_wp, -0.405648_wp, -0.237302_wp, &
                      0.139178_wp, 0.131263_wp, 0.139578_wp, -0.038152_wp, 0.142137_wp, &
                      0.115989_wp, 0.128096_wp, 0.20004_wp, 0.189686_wp]
         case ("mbis")
            charges = [ &
                      -0.471257_wp, 0.037183_wp, -0.463274_wp, -0.202543_wp, -0.260933_wp, &
                      0.147506_wp, 0.145782_wp, 0.151222_wp, 0.122358_wp, 0.142982_wp, &
                      0.151589_wp, 0.145029_wp, 0.176946_wp, 0.177419_wp]
         case ("chelpg")
            charges = [ &
                      -0.598471_wp, 0.642991_wp, -0.571719_wp, -0.418021_wp, -0.235999_wp, &
                      0.144221_wp, 0.135878_wp, 0.144426_wp, -0.039577_wp, 0.146768_wp, &
                      0.12058_wp, 0.132468_wp, 0.203528_wp, 0.192929_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (9)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070552_wp, -0.034403_wp, -0.036376_wp, -0.036276_wp, -0.034539_wp, &
                      -0.033767_wp, -0.036835_wp, 0.013076_wp, -0.139337_wp, 0.022113_wp, &
                      0.023946_wp, 0.022108_wp, 0.01974_wp, 0.01974_wp, 0.019447_wp, &
                      0.019464_wp, 0.019524_wp, 0.019525_wp, 0.021541_wp, 0.021552_wp, &
                      0.021798_wp, 0.021809_wp, 0.028435_wp, 0.028462_wp, 0.039939_wp, &
                      0.039913_wp]
         case ("resp")
            charges = [ &
                      -0.337726_wp, 0.197662_wp, 0.005699_wp, -0.012476_wp, 0.032888_wp, &
                      0.031682_wp, 0.11965_wp, -0.270947_wp, -0.182407_wp, 0.072375_wp, &
                      0.080315_wp, 0.073986_wp, -0.031128_wp, -0.033093_wp, -7.6e-05_wp, &
                      -0.000518_wp, -0.009383_wp, -0.006955_wp, -0.015411_wp, -0.013192_wp, &
                      -0.005182_wp, -0.006115_wp, 0.011158_wp, 0.010087_wp, 0.144253_wp, &
                      0.144856_wp]
         case ("mbis")
            charges = [ &
                      -0.421193_wp, -0.160801_wp, -0.211721_wp, -0.214728_wp, -0.211416_wp, &
                      -0.212071_wp, -0.221318_wp, -0.188086_wp, -0.203351_wp, 0.128147_wp, &
                      0.130376_wp, 0.128086_wp, 0.100365_wp, 0.100411_wp, 0.107084_wp, &
                      0.107074_wp, 0.106915_wp, 0.106903_wp, 0.109845_wp, 0.109855_wp, &
                      0.110686_wp, 0.110666_wp, 0.136005_wp, 0.135876_wp, 0.158223_wp, &
                      0.158216_wp]
         case ("chelpg")
            charges = [ &
                      -0.354198_wp, 0.211245_wp, 0.000638_wp, -0.013338_wp, 0.039079_wp, &
                      0.025948_wp, 0.129476_wp, -0.280482_wp, -0.18132_wp, 0.076106_wp, &
                      0.084032_wp, 0.077699_wp, -0.033675_wp, -0.035647_wp, 0.000631_wp, &
                      0.000174_wp, -0.009407_wp, -0.006971_wp, -0.016867_wp, -0.014633_wp, &
                      -0.004534_wp, -0.005466_wp, 0.009246_wp, 0.008161_wp, 0.146755_wp, &
                      0.147349_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073014_wp, -0.035171_wp, -0.037046_wp, -0.036638_wp, -0.034583_wp, &
                      -0.031763_wp, -0.036124_wp, 0.019601_wp, -0.171558_wp, 0.0222_wp, &
                      0.023639_wp, 0.022201_wp, 0.020159_wp, 0.02016_wp, 0.0198_wp, &
                      0.019827_wp, 0.020249_wp, 0.020264_wp, 0.022196_wp, 0.022214_wp, &
                      0.024973_wp, 0.024983_wp, 0.029119_wp, 0.029164_wp, 0.047625_wp, &
                      0.04757_wp]
         case ("resp")
            charges = [ &
                      -0.356018_wp, 0.219979_wp, 0.012061_wp, -0.014619_wp, 0.04454_wp, &
                      0.038965_wp, 0.114352_wp, -0.221452_wp, -0.238758_wp, 0.074686_wp, &
                      0.082111_wp, 0.076315_wp, -0.03668_wp, -0.038692_wp, -0.003354_wp, &
                      -0.003837_wp, -0.010727_wp, -0.008234_wp, -0.018623_wp, -0.016305_wp, &
                      -0.0042_wp, -0.005089_wp, 0.010203_wp, 0.009102_wp, 0.146866_wp, &
                      0.147409_wp]
         case ("mbis")
            charges = [ &
                      -0.43311_wp, -0.160763_wp, -0.214461_wp, -0.217564_wp, -0.214426_wp, &
                      -0.211653_wp, -0.225496_wp, -0.175084_wp, -0.245863_wp, 0.130763_wp, &
                      0.13229_wp, 0.130714_wp, 0.101533_wp, 0.101578_wp, 0.108394_wp, &
                      0.108402_wp, 0.108701_wp, 0.108717_wp, 0.111422_wp, 0.111444_wp, &
                      0.11554_wp, 0.115517_wp, 0.13756_wp, 0.137456_wp, 0.16925_wp, &
                      0.16919_wp]
         case ("chelpg")
            charges = [ &
                      -0.372882_wp, 0.233632_wp, 0.008004_wp, -0.016684_wp, 0.05188_wp, &
                      0.033394_wp, 0.123586_wp, -0.230513_wp, -0.237718_wp, 0.078488_wp, &
                      0.085942_wp, 0.080099_wp, -0.03926_wp, -0.041278_wp, -0.002878_wp, &
                      -0.003377_wp, -0.010585_wp, -0.008083_wp, -0.020345_wp, -0.018011_wp, &
                      -0.003615_wp, -0.004503_wp, 0.008381_wp, 0.007268_wp, 0.149262_wp, &
                      0.149795_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073709_wp, -0.035368_wp, -0.037231_wp, -0.036666_wp, -0.034559_wp, &
                      -0.030989_wp, -0.035851_wp, 0.0219_wp, -0.182043_wp, 0.022169_wp, &
                      0.023466_wp, 0.022171_wp, 0.020279_wp, 0.02028_wp, 0.019832_wp, &
                      0.019862_wp, 0.020519_wp, 0.020538_wp, 0.022289_wp, 0.02231_wp, &
                      0.02613_wp, 0.026136_wp, 0.029168_wp, 0.029221_wp, 0.050129_wp, &
                      0.050064_wp]
         case ("resp")
            charges = [ &
                      -0.361117_wp, 0.226897_wp, 0.013744_wp, -0.015493_wp, 0.048585_wp, &
                      0.041641_wp, 0.10981_wp, -0.201192_wp, -0.257973_wp, 0.075226_wp, &
                      0.082422_wp, 0.07686_wp, -0.038395_wp, -0.040419_wp, -0.004437_wp, &
                      -0.00493_wp, -0.010947_wp, -0.008438_wp, -0.019786_wp, -0.017439_wp, &
                      -0.003639_wp, -0.00452_wp, 0.010124_wp, 0.009018_wp, 0.146937_wp, &
                      0.147459_wp]
         case ("mbis")
            charges = [ &
                      -0.436268_wp, -0.160709_wp, -0.215143_wp, -0.21825_wp, -0.215114_wp, &
                      -0.2114_wp, -0.226581_wp, -0.170209_wp, -0.259774_wp, 0.131386_wp, &
                      0.132657_wp, 0.13134_wp, 0.101853_wp, 0.101898_wp, 0.108628_wp, &
                      0.108643_wp, 0.109277_wp, 0.109301_wp, 0.111662_wp, 0.111689_wp, &
                      0.117284_wp, 0.117255_wp, 0.137707_wp, 0.137615_wp, 0.172689_wp, &
                      0.172612_wp]
         case ("chelpg")
            charges = [ &
                      -0.378103_wp, 0.240563_wp, 0.01002_wp, -0.017981_wp, 0.056299_wp, &
                      0.036177_wp, 0.118729_wp, -0.210015_wp, -0.256957_wp, 0.07905_wp, &
                      0.086289_wp, 0.080666_wp, -0.040981_wp, -0.043011_wp, -0.004034_wp, &
                      -0.004544_wp, -0.010745_wp, -0.008226_wp, -0.021593_wp, -0.01923_wp, &
                      -0.003082_wp, -0.003962_wp, 0.008358_wp, 0.007241_wp, 0.149279_wp, &
                      0.149792_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (10)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069609_wp, -0.032156_wp, -0.033762_wp, -0.036776_wp, 0.013112_wp, &
                      -0.138951_wp, 0.023088_wp, 0.025737_wp, 0.023095_wp, 0.022627_wp, &
                      0.022634_wp, 0.021953_wp, 0.021945_wp, 0.02856_wp, 0.028559_wp, &
                      0.03996_wp, 0.039996_wp]
         case ("resp")
            charges = [ &
                      -0.307436_wp, 0.184427_wp, -0.01282_wp, 0.075254_wp, -0.208184_wp, &
                      -0.188699_wp, 0.066964_wp, 0.077386_wp, 0.068131_wp, -0.029005_wp, &
                      -0.030318_wp, 0.006773_wp, 0.007744_wp, 0.015003_wp, 0.01691_wp, &
                      0.128996_wp, 0.128872_wp]
         case ("mbis")
            charges = [ &
                      -0.420724_wp, -0.163489_wp, -0.213223_wp, -0.225297_wp, -0.187117_wp, &
                      -0.202958_wp, 0.129186_wp, 0.132938_wp, 0.129132_wp, 0.104917_wp, &
                      0.104986_wp, 0.111353_wp, 0.111284_wp, 0.136338_wp, 0.136353_wp, &
                      0.158171_wp, 0.158165_wp]
         case ("chelpg")
            charges = [ &
                      -0.321666_wp, 0.198345_wp, -0.022933_wp, 0.084703_wp, -0.216435_wp, &
                      -0.187883_wp, 0.070257_wp, 0.080425_wp, 0.071418_wp, -0.03176_wp, &
                      -0.033085_wp, 0.008404_wp, 0.009369_wp, 0.01325_wp, 0.015181_wp, &
                      0.13126_wp, 0.131148_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071757_wp, -0.032141_wp, -0.031494_wp, -0.035945_wp, 0.020094_wp, &
                      -0.173239_wp, 0.02346_wp, 0.02515_wp, 0.023475_wp, 0.023345_wp, &
                      0.023329_wp, 0.025466_wp, 0.025446_wp, 0.029325_wp, 0.029302_wp, &
                      0.048082_wp, 0.048115_wp]
         case ("resp")
            charges = [ &
                      -0.330349_wp, 0.205456_wp, -0.006321_wp, 0.065204_wp, -0.154084_wp, &
                      -0.248772_wp, 0.070846_wp, 0.079536_wp, 0.072045_wp, -0.03303_wp, &
                      -0.034401_wp, 0.009005_wp, 0.009895_wp, 0.015124_wp, 0.017022_wp, &
                      0.131424_wp, 0.131401_wp]
         case ("mbis")
            charges = [ &
                      -0.432805_wp, -0.163525_wp, -0.212814_wp, -0.22972_wp, -0.173136_wp, &
                      -0.248208_wp, 0.1322_wp, 0.13433_wp, 0.132163_wp, 0.106295_wp, &
                      0.10632_wp, 0.116732_wp, 0.116637_wp, 0.138014_wp, 0.137989_wp, &
                      0.169775_wp, 0.169766_wp]
         case ("chelpg")
            charges = [ &
                      -0.344945_wp, 0.219462_wp, -0.015747_wp, 0.073667_wp, -0.161627_wp, &
                      -0.248028_wp, 0.074208_wp, 0.082683_wp, 0.075402_wp, -0.035814_wp, &
                      -0.037196_wp, 0.010492_wp, 0.011377_wp, 0.013556_wp, 0.015477_wp, &
                      0.133522_wp, 0.133512_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072195_wp, -0.032118_wp, -0.030835_wp, -0.035702_wp, 0.021971_wp, &
                      -0.181749_wp, 0.023557_wp, 0.024864_wp, 0.023574_wp, 0.023425_wp, &
                      0.023402_wp, 0.026434_wp, 0.026411_wp, 0.029381_wp, 0.029351_wp, &
                      0.050105_wp, 0.050137_wp]
         case ("resp")
            charges = [ &
                      -0.335664_wp, 0.210481_wp, -0.004139_wp, 0.060542_wp, -0.137571_wp, &
                      -0.26435_wp, 0.071808_wp, 0.079743_wp, 0.073012_wp, -0.034101_wp, &
                      -0.035488_wp, 0.009676_wp, 0.010543_wp, 0.015339_wp, 0.017229_wp, &
                      0.13147_wp, 0.131471_wp]
         case ("mbis")
            charges = [ &
                      -0.435211_wp, -0.163443_wp, -0.212606_wp, -0.230588_wp, -0.169151_wp, &
                      -0.2595_wp, 0.132857_wp, 0.134317_wp, 0.132821_wp, 0.10644_wp, &
                      0.106453_wp, 0.118195_wp, 0.118093_wp, 0.138149_wp, 0.138113_wp, &
                      0.172543_wp, 0.172532_wp]
         case ("chelpg")
            charges = [ &
                      -0.350366_wp, 0.224468_wp, -0.013255_wp, 0.06854_wp, -0.144785_wp, &
                      -0.263638_wp, 0.075189_wp, 0.082925_wp, 0.076387_wp, -0.036883_wp, &
                      -0.038281_wp, 0.011104_wp, 0.011965_wp, 0.013861_wp, 0.015773_wp, &
                      0.133491_wp, 0.133504_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (11)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.065118_wp, -0.033068_wp, 0.013776_wp, -0.139047_wp, 0.031079_wp, &
                      0.026152_wp, 0.026136_wp, 0.029873_wp, 0.029885_wp, 0.040162_wp, &
                      0.040175_wp]
         case ("resp")
            charges = [ &
                      -0.344091_wp, 0.295789_wp, -0.307483_wp, -0.18484_wp, 0.083173_wp, &
                      0.084859_wp, 0.085553_wp, -0.012959_wp, -0.012881_wp, 0.156886_wp, &
                      0.155996_wp]
         case ("mbis")
            charges = [ &
                      -0.422181_wp, -0.173817_wp, -0.188213_wp, -0.203884_wp, 0.141023_wp, &
                      0.133965_wp, 0.133921_wp, 0.131062_wp, 0.131038_wp, 0.158557_wp, &
                      0.158534_wp]
         case ("chelpg")
            charges = [ &
                      -0.360807_wp, 0.30991_wp, -0.318032_wp, -0.183799_wp, 0.08685_wp, &
                      0.088757_wp, 0.089445_wp, -0.015327_wp, -0.01524_wp, 0.159562_wp, &
                      0.158683_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.064007_wp, -0.032015_wp, 0.021252_wp, -0.175385_wp, 0.032775_wp, &
                      0.029235_wp, 0.029226_wp, 0.030645_wp, 0.030659_wp, 0.048806_wp, &
                      0.048817_wp]
         case ("resp")
            charges = [ &
                      -0.360096_wp, 0.307333_wp, -0.261548_wp, -0.248081_wp, 0.088211_wp, &
                      0.090702_wp, 0.091422_wp, -0.016255_wp, -0.016103_wp, 0.162648_wp, &
                      0.161768_wp]
         case ("mbis")
            charges = [ &
                      -0.429656_wp, -0.175309_wp, -0.173774_wp, -0.251762_wp, 0.144546_wp, &
                      0.139741_wp, 0.139711_wp, 0.132362_wp, 0.13234_wp, 0.170916_wp, &
                      0.170891_wp]
         case ("chelpg")
            charges = [ &
                      -0.377005_wp, 0.321586_wp, -0.27205_wp, -0.247053_wp, 0.091939_wp, &
                      0.094642_wp, 0.095355_wp, -0.018655_wp, -0.018494_wp, 0.165302_wp, &
                      0.164433_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.063688_wp, -0.031801_wp, 0.02274_wp, -0.182119_wp, 0.032999_wp, &
                      0.029848_wp, 0.029841_wp, 0.030676_wp, 0.03069_wp, 0.050405_wp, &
                      0.050416_wp]
         case ("resp")
            charges = [ &
                      -0.362069_wp, 0.307608_wp, -0.250665_wp, -0.2603_wp, 0.088909_wp, &
                      0.09169_wp, 0.092414_wp, -0.01671_wp, -0.016545_wp, 0.163273_wp, &
                      0.162395_wp]
         case ("mbis")
            charges = [ &
                      -0.430565_wp, -0.17552_wp, -0.17071_wp, -0.260674_wp, 0.144944_wp, &
                      0.1408_wp, 0.140773_wp, 0.13239_wp, 0.132369_wp, 0.173112_wp, &
                      0.173087_wp]
         case ("chelpg")
            charges = [ &
                      -0.37902_wp, 0.321859_wp, -0.261132_wp, -0.25928_wp, 0.09265_wp, &
                      0.09564_wp, 0.096357_wp, -0.01911_wp, -0.018936_wp, 0.165919_wp, &
                      0.165052_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (12)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068594_wp, -0.032034_wp, -0.043137_wp, 0.041577_wp, -0.250832_wp, &
                      0.024925_wp, 0.023579_wp, 0.022473_wp, 0.020388_wp, 0.02224_wp, &
                      0.018424_wp, 0.024384_wp, 0.031492_wp, 0.018155_wp, 0.146963_wp]
         case ("resp")
            charges = [ &
                      -0.272235_wp, 0.198428_wp, -0.130558_wp, 0.203361_wp, -0.613584_wp, &
                      0.066354_wp, 0.063524_wp, 0.058983_wp, -0.031319_wp, -0.038107_wp, &
                      0.032565_wp, 0.046025_wp, 0.057645_wp, -0.013954_wp, 0.372874_wp]
         case ("mbis")
            charges = [ &
                      -0.422158_wp, -0.135404_wp, -0.323077_wp, 0.091276_wp, -0.612793_wp, &
                      0.131288_wp, 0.130845_wp, 0.128874_wp, 0.098083_wp, 0.101065_wp, &
                      0.11342_wp, 0.135047_wp, 0.093852_wp, 0.056694_wp, 0.412991_wp]
         case ("chelpg")
            charges = [ &
                      -0.287987_wp, 0.214363_wp, -0.145192_wp, 0.21343_wp, -0.61569_wp, &
                      0.069581_wp, 0.06717_wp, 0.062654_wp, -0.034024_wp, -0.041004_wp, &
                      0.035158_wp, 0.04895_wp, 0.05557_wp, -0.016393_wp, 0.373413_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070825_wp, -0.031353_wp, -0.042547_wp, 0.04019_wp, -0.279093_wp, &
                      0.025216_wp, 0.023532_wp, 0.023443_wp, 0.023197_wp, 0.023573_wp, &
                      0.025195_wp, 0.021783_wp, 0.030729_wp, 0.024347_wp, 0.162613_wp]
         case ("resp")
            charges = [ &
                      -0.293436_wp, 0.212628_wp, -0.118239_wp, 0.235267_wp, -0.707278_wp, &
                      0.070237_wp, 0.066662_wp, 0.062728_wp, -0.031892_wp, -0.040215_wp, &
                      0.041152_wp, 0.035139_wp, 0.048833_wp, -0.007131_wp, 0.425546_wp]
         case ("mbis")
            charges = [ &
                      -0.436008_wp, -0.134893_wp, -0.325515_wp, 0.099259_wp, -0.685521_wp, &
                      0.134385_wp, 0.133321_wp, 0.133042_wp, 0.102757_wp, 0.103256_wp, &
                      0.126071_wp, 0.130187_wp, 0.092908_wp, 0.067625_wp, 0.459127_wp]
         case ("chelpg")
            charges = [ &
                      -0.309409_wp, 0.228676_wp, -0.132739_wp, 0.245373_wp, -0.709381_wp, &
                      0.073521_wp, 0.070361_wp, 0.066448_wp, -0.034631_wp, -0.043139_wp, &
                      0.043702_wp, 0.038009_wp, 0.046731_wp, -0.009594_wp, 0.426072_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070966_wp, -0.031282_wp, -0.042482_wp, 0.040097_wp, -0.281115_wp, &
                      0.025218_wp, 0.023481_wp, 0.023537_wp, 0.023463_wp, 0.023624_wp, &
                      0.025774_wp, 0.021502_wp, 0.030562_wp, 0.024874_wp, 0.163712_wp]
         case ("resp")
            charges = [ &
                      -0.294901_wp, 0.213705_wp, -0.11724_wp, 0.237519_wp, -0.714053_wp, &
                      0.070468_wp, 0.066788_wp, 0.063041_wp, -0.031817_wp, -0.040467_wp, &
                      0.041943_wp, 0.034139_wp, 0.047984_wp, -0.00643_wp, 0.42932_wp]
         case ("mbis")
            charges = [ &
                      -0.436885_wp, -0.134845_wp, -0.325595_wp, 0.099885_wp, -0.690699_wp, &
                      0.134555_wp, 0.133394_wp, 0.133366_wp, 0.1032_wp, 0.103321_wp, &
                      0.127135_wp, 0.12964_wp, 0.092627_wp, 0.068559_wp, 0.462343_wp]
         case ("chelpg")
            charges = [ &
                      -0.310891_wp, 0.229754_wp, -0.131724_wp, 0.247618_wp, -0.716159_wp, &
                      0.073757_wp, 0.070491_wp, 0.066764_wp, -0.034555_wp, -0.043392_wp, &
                      0.044492_wp, 0.037003_wp, 0.04588_wp, -0.008889_wp, 0.42985_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (13)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.02447_wp, -0.036611_wp, -0.033255_wp, -0.034234_wp, -0.033567_wp, &
                      -0.069404_wp, 0.021399_wp, 0.021425_wp, 0.028415_wp, 0.028396_wp, &
                      0.039314_wp, 0.039332_wp, -0.148139_wp, 0.021128_wp, 0.021147_wp, &
                      0.020117_wp, 0.020103_wp, 0.022603_wp, 0.02261_wp, 0.024762_wp]
         case ("resp")
            charges = [ &
                      -0.211262_wp, 0.173302_wp, -0.074397_wp, -0.021071_wp, 0.235198_wp, &
                      -0.364719_wp, 0.009532_wp, 0.011203_wp, 0.00036_wp, 0.00054_wp, &
                      0.130321_wp, 0.129488_wp, -0.202133_wp, 0.00311_wp, 0.004013_wp, &
                      -0.036433_wp, -0.036114_wp, 0.082324_wp, 0.081311_wp, 0.085429_wp]
         case ("mbis")
            charges = [ &
                      -0.117596_wp, -0.230586_wp, -0.217051_wp, -0.212741_wp, -0.161338_wp, &
                      -0.421394_wp, 0.110272_wp, 0.110317_wp, 0.134755_wp, 0.134751_wp, &
                      0.145403_wp, 0.145444_wp, -0.23292_wp, 0.110152_wp, 0.110176_wp, &
                      0.101382_wp, 0.101364_wp, 0.128919_wp, 0.128943_wp, 0.131756_wp]
         case ("chelpg")
            charges = [ &
                      -0.224593_wp, 0.190162_wp, -0.085696_wp, -0.023744_wp, 0.251093_wp, &
                      -0.384415_wp, 0.011543_wp, 0.013207_wp, -0.00258_wp, -0.002409_wp, &
                      0.133596_wp, 0.13275_wp, -0.200734_wp, 0.00344_wp, 0.004344_wp, &
                      -0.039235_wp, -0.038922_wp, 0.086695_wp, 0.085676_wp, 0.089821_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.029457_wp, -0.036129_wp, -0.031202_wp, -0.034304_wp, -0.033993_wp, &
                      -0.071938_wp, 0.024798_wp, 0.024806_wp, 0.028931_wp, 0.028911_wp, &
                      0.047077_wp, 0.047091_wp, -0.178632_wp, 0.021868_wp, 0.021876_wp, &
                      0.020908_wp, 0.020901_wp, 0.02259_wp, 0.022589_wp, 0.024405_wp]
         case ("resp")
            charges = [ &
                      -0.174547_wp, 0.173303_wp, -0.067424_wp, -0.020471_wp, 0.252235_wp, &
                      -0.381587_wp, 0.011305_wp, 0.013001_wp, -0.001561_wp, -0.00143_wp, &
                      0.136426_wp, 0.135572_wp, -0.256113_wp, 0.002756_wp, 0.003619_wp, &
                      -0.040046_wp, -0.039731_wp, 0.084237_wp, 0.083211_wp, 0.087243_wp]
         case ("mbis")
            charges = [ &
                      -0.109474_wp, -0.234177_wp, -0.216492_wp, -0.216277_wp, -0.161216_wp, &
                      -0.433657_wp, 0.115487_wp, 0.115497_wp, 0.135986_wp, 0.135981_wp, &
                      0.157421_wp, 0.157452_wp, -0.273022_wp, 0.111962_wp, 0.111962_wp, &
                      0.103016_wp, 0.10301_wp, 0.131431_wp, 0.131441_wp, 0.133676_wp]
         case ("chelpg")
            charges = [ &
                      -0.187679_wp, 0.189986_wp, -0.078421_wp, -0.023288_wp, 0.268306_wp, &
                      -0.40139_wp, 0.013243_wp, 0.014932_wp, -0.004484_wp, -0.004362_wp, &
                      0.139643_wp, 0.138776_wp, -0.254733_wp, 0.003098_wp, 0.003962_wp, &
                      -0.042892_wp, -0.042583_wp, 0.08863_wp, 0.087599_wp, 0.091656_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.030856_wp, -0.035976_wp, -0.030555_wp, -0.034286_wp, -0.034032_wp, &
                      -0.07254_wp, 0.025792_wp, 0.025793_wp, 0.028903_wp, 0.028883_wp, &
                      0.049087_wp, 0.0491_wp, -0.18646_wp, 0.021956_wp, 0.02196_wp, &
                      0.02114_wp, 0.021137_wp, 0.02251_wp, 0.022507_wp, 0.024232_wp]
         case ("resp")
            charges = [ &
                      -0.163099_wp, 0.171431_wp, -0.064627_wp, -0.020437_wp, 0.256201_wp, &
                      -0.38535_wp, 0.011868_wp, 0.013566_wp, -0.001984_wp, -0.001865_wp, &
                      0.137682_wp, 0.136822_wp, -0.270383_wp, 0.002565_wp, 0.003416_wp, &
                      -0.040814_wp, -0.040498_wp, 0.084521_wp, 0.08349_wp, 0.087496_wp]
         case ("mbis")
            charges = [ &
                      -0.106989_wp, -0.23483_wp, -0.216234_wp, -0.21693_wp, -0.161125_wp, &
                      -0.436331_wp, 0.116993_wp, 0.116991_wp, 0.135976_wp, 0.13597_wp, &
                      0.160438_wp, 0.160466_wp, -0.283381_wp, 0.112198_wp, 0.112191_wp, &
                      0.103454_wp, 0.103454_wp, 0.13187_wp, 0.131876_wp, 0.133952_wp]
         case ("chelpg")
            charges = [ &
                      -0.176113_wp, 0.187991_wp, -0.075488_wp, -0.023329_wp, 0.272303_wp, &
                      -0.40515_wp, 0.013777_wp, 0.015469_wp, -0.004888_wp, -0.004777_wp, &
                      0.140869_wp, 0.139996_wp, -0.269021_wp, 0.00292_wp, 0.003771_wp, &
                      -0.043669_wp, -0.043357_wp, 0.088913_wp, 0.087877_wp, 0.091905_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (14)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.033385_wp, -0.031783_wp, -0.036627_wp, 0.022207_wp, 0.022211_wp, &
                      0.021439_wp, 0.021446_wp, 0.028483_wp, 0.028471_wp, -0.069053_wp, &
                      0.022866_wp, 0.02286_wp, 0.025758_wp, 0.024558_wp, 0.039355_wp, &
                      0.039374_wp, -0.148184_wp]
         case ("resp")
            charges = [ &
                      -0.048577_wp, 0.188433_wp, 0.059161_wp, -0.02711_wp, -0.027013_wp, &
                      0.015815_wp, 0.015204_wp, 0.019947_wp, 0.019679_wp, -0.313894_wp, &
                      0.070636_wp, 0.070477_wp, 0.078665_wp, -0.143477_wp, 0.113337_wp, &
                      0.114402_wp, -0.205686_wp]
         case ("mbis")
            charges = [ &
                      -0.214328_wp, -0.162117_wp, -0.235903_wp, 0.104621_wp, 0.104633_wp, &
                      0.110577_wp, 0.110596_wp, 0.135217_wp, 0.135295_wp, -0.421626_wp, &
                      0.129196_wp, 0.129177_wp, 0.133311_wp, -0.116435_wp, 0.145334_wp, &
                      0.145357_wp, -0.232908_wp]
         case ("chelpg")
            charges = [ &
                      -0.061399_wp, 0.20377_wp, 0.06883_wp, -0.030012_wp, -0.029899_wp, &
                      0.018122_wp, 0.017518_wp, 0.018255_wp, 0.017986_wp, -0.328917_wp, &
                      0.074096_wp, 0.073958_wp, 0.081786_wp, -0.151006_wp, 0.115408_wp, &
                      0.116462_wp, -0.204958_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.031237_wp, -0.031891_wp, -0.036103_wp, 0.0229_wp, 0.022901_wp, &
                      0.024948_wp, 0.024948_wp, 0.02905_wp, 0.029024_wp, -0.071236_wp, &
                      0.023268_wp, 0.023265_wp, 0.025305_wp, 0.029669_wp, 0.04727_wp, &
                      0.047289_wp, -0.179372_wp]
         case ("resp")
            charges = [ &
                      -0.039119_wp, 0.20443_wp, 0.052208_wp, -0.029971_wp, -0.029848_wp, &
                      0.017895_wp, 0.017297_wp, 0.019297_wp, 0.019009_wp, -0.336562_wp, &
                      0.074601_wp, 0.074461_wp, 0.081446_wp, -0.103174_wp, 0.118904_wp, &
                      0.120004_wp, -0.260879_wp]
         case ("mbis")
            charges = [ &
                      -0.213713_wp, -0.162607_wp, -0.239687_wp, 0.106036_wp, 0.106042_wp, &
                      0.115976_wp, 0.11598_wp, 0.136539_wp, 0.136596_wp, -0.43379_wp, &
                      0.132259_wp, 0.132246_wp, 0.134947_wp, -0.108105_wp, 0.157583_wp, &
                      0.157605_wp, -0.273912_wp]
         case ("chelpg")
            charges = [ &
                      -0.05114_wp, 0.21975_wp, 0.06083_wp, -0.032874_wp, -0.032735_wp, &
                      0.020034_wp, 0.019445_wp, 0.017782_wp, 0.017494_wp, -0.351941_wp, &
                      0.078131_wp, 0.078012_wp, 0.084677_wp, -0.109812_wp, 0.120756_wp, &
                      0.121845_wp, -0.260254_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.030638_wp, -0.031894_wp, -0.035954_wp, 0.022967_wp, 0.022968_wp, &
                      0.025864_wp, 0.025862_wp, 0.029038_wp, 0.029009_wp, -0.07166_wp, &
                      0.023366_wp, 0.023365_wp, 0.025067_wp, 0.030955_wp, 0.049114_wp, &
                      0.049133_wp, -0.186564_wp]
         case ("resp")
            charges = [ &
                      -0.036103_wp, 0.207705_wp, 0.049155_wp, -0.030673_wp, -0.030545_wp, &
                      0.018467_wp, 0.017874_wp, 0.01917_wp, 0.018879_wp, -0.341353_wp, &
                      0.075481_wp, 0.075345_wp, 0.081776_wp, -0.092224_wp, 0.119956_wp, &
                      0.121063_wp, -0.273972_wp]
         case ("mbis")
            charges = [ &
                      -0.213465_wp, -0.162618_wp, -0.240305_wp, 0.106176_wp, 0.106181_wp, &
                      0.117367_wp, 0.117368_wp, 0.13655_wp, 0.136602_wp, -0.436099_wp, &
                      0.132893_wp, 0.132882_wp, 0.134996_wp, -0.105826_wp, 0.160351_wp, &
                      0.160373_wp, -0.283428_wp]
         case ("chelpg")
            charges = [ &
                      -0.047818_wp, 0.222972_wp, 0.057345_wp, -0.033566_wp, -0.033421_wp, &
                      0.020547_wp, 0.019962_wp, 0.017734_wp, 0.017442_wp, -0.356821_wp, &
                      0.079027_wp, 0.078914_wp, 0.08504_wp, -0.098501_wp, 0.12172_wp, &
                      0.122817_wp, -0.273392_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (15)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.033676_wp, -0.065346_wp, 0.024173_wp, 0.026047_wp, 0.026048_wp, &
                      0.030931_wp, 0.02984_wp, 0.02984_wp, 0.03949_wp, 0.039491_wp, &
                      -0.146835_wp]
         case ("resp")
            charges = [ &
                      0.257423_wp, -0.352187_wp, -0.235507_wp, 0.087243_wp, 0.087175_wp, &
                      0.087388_wp, -0.004087_wp, -0.004781_wp, 0.139414_wp, 0.139637_wp, &
                      -0.201717_wp]
         case ("mbis")
            charges = [ &
                      -0.185225_wp, -0.42309_wp, -0.117725_wp, 0.133443_wp, 0.133441_wp, &
                      0.141299_wp, 0.13022_wp, 0.130224_wp, 0.145598_wp, 0.145596_wp, &
                      -0.233776_wp]
         case ("chelpg")
            charges = [ &
                      0.271039_wp, -0.368922_wp, -0.245762_wp, 0.091138_wp, 0.091088_wp, &
                      0.091125_wp, -0.006323_wp, -0.007031_wp, 0.142025_wp, 0.142242_wp, &
                      -0.200619_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.032957_wp, -0.064288_wp, 0.029549_wp, 0.02918_wp, 0.02918_wp, &
                      0.032523_wp, 0.030372_wp, 0.030371_wp, 0.047792_wp, 0.047793_wp, &
                      -0.179511_wp]
         case ("resp")
            charges = [ &
                      0.272985_wp, -0.365173_wp, -0.202972_wp, 0.092797_wp, 0.092752_wp, &
                      0.091047_wp, -0.008409_wp, -0.009143_wp, 0.147744_wp, 0.148007_wp, &
                      -0.259635_wp]
         case ("mbis")
            charges = [ &
                      -0.185951_wp, -0.430235_wp, -0.109364_wp, 0.139311_wp, 0.13931_wp, &
                      0.14456_wp, 0.131094_wp, 0.131098_wp, 0.158483_wp, 0.158481_wp, &
                      -0.276783_wp]
         case ("chelpg")
            charges = [ &
                      0.286772_wp, -0.382099_wp, -0.213174_wp, 0.096731_wp, 0.096704_wp, &
                      0.09483_wp, -0.010693_wp, -0.01144_wp, 0.150331_wp, 0.150588_wp, &
                      -0.25855_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.032826_wp, -0.06401_wp, 0.030547_wp, 0.029773_wp, 0.029773_wp, &
                      0.032722_wp, 0.03035_wp, 0.03035_wp, 0.04923_wp, 0.049231_wp, &
                      -0.185135_wp]
         case ("resp")
            charges = [ &
                      0.274458_wp, -0.366428_wp, -0.196078_wp, 0.093676_wp, 0.093635_wp, &
                      0.091397_wp, -0.009135_wp, -0.009875_wp, 0.148979_wp, 0.149248_wp, &
                      -0.269878_wp]
         case ("mbis")
            charges = [ &
                      -0.185988_wp, -0.43108_wp, -0.107661_wp, 0.140342_wp, 0.14034_wp, &
                      0.144903_wp, 0.131032_wp, 0.131036_wp, 0.160653_wp, 0.160652_wp, &
                      -0.284224_wp]
         case ("chelpg")
            charges = [ &
                      0.288264_wp, -0.383386_wp, -0.206255_wp, 0.097617_wp, 0.097593_wp, &
                      0.09519_wp, -0.011425_wp, -0.012177_wp, 0.151557_wp, 0.15182_wp, &
                      -0.268798_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (16)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.250903_wp, 0.0415_wp, -0.043042_wp, -0.034286_wp, -0.035172_wp, &
                      -0.035996_wp, -0.036007_wp, -0.036121_wp, -0.036343_wp, -0.03411_wp, &
                      -0.070246_wp, 0.146826_wp, 0.018046_wp, 0.031353_wp, 0.024242_wp, &
                      0.018319_wp, 0.019256_wp, 0.021092_wp, 0.019901_wp, 0.018802_wp, &
                      0.018462_wp, 0.018883_wp, 0.018698_wp, 0.018378_wp, 0.018315_wp, &
                      0.018493_wp, 0.01841_wp, 0.018276_wp, 0.01915_wp, 0.019245_wp, &
                      0.021589_wp, 0.021509_wp, 0.023499_wp]
         case ("resp")
            charges = [ &
                      -0.632256_wp, 0.245343_wp, -0.067646_wp, -0.002025_wp, 0.081568_wp, &
                      0.003546_wp, 0.060477_wp, -0.044205_wp, -0.001042_wp, 0.19541_wp, &
                      -0.335701_wp, 0.375477_wp, -0.022806_wp, 0.050508_wp, 0.028549_wp, &
                      0.012739_wp, -0.006216_wp, -0.016257_wp, -0.020054_wp, -0.023876_wp, &
                      -0.007532_wp, -0.002431_wp, -0.018405_wp, -0.018311_wp, -0.000195_wp, &
                      -0.001795_wp, 0.001799_wp, 0.001276_wp, -0.029892_wp, -0.029802_wp, &
                      0.072359_wp, 0.072283_wp, 0.079113_wp]
         case ("mbis")
            charges = [ &
                      -0.613547_wp, 0.095104_wp, -0.32322_wp, -0.182975_wp, -0.211924_wp, &
                      -0.212122_wp, -0.207917_wp, -0.214153_wp, -0.211358_wp, -0.160202_wp, &
                      -0.421643_wp, 0.412935_wp, 0.056342_wp, 0.09345_wp, 0.134769_wp, &
                      0.113062_wp, 0.102748_wp, 0.105737_wp, 0.108322_wp, 0.10636_wp, &
                      0.105531_wp, 0.106263_wp, 0.105546_wp, 0.104926_wp, 0.105436_wp, &
                      0.105742_wp, 0.105924_wp, 0.105642_wp, 0.099908_wp, 0.100075_wp, &
                      0.1277_wp, 0.127531_wp, 0.130026_wp]
         case ("chelpg")
            charges = [ &
                      -0.634644_wp, 0.255702_wp, -0.079005_wp, 0.001292_wp, 0.089217_wp, &
                      -0.008049_wp, 0.074983_wp, -0.051664_wp, -0.005232_wp, 0.210012_wp, &
                      -0.353511_wp, 0.375975_wp, -0.025222_wp, 0.048464_wp, 0.030892_wp, &
                      0.014748_wp, -0.006572_wp, -0.01686_wp, -0.02186_wp, -0.025607_wp, &
                      -0.005948_wp, -0.000739_wp, -0.021158_wp, -0.021063_wp, 0.001182_wp, &
                      -0.000449_wp, 0.002414_wp, 0.001872_wp, -0.032505_wp, -0.032414_wp, &
                      0.076364_wp, 0.07629_wp, 0.083094_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.27657_wp, 0.040178_wp, -0.04262_wp, -0.033828_wp, -0.035701_wp, &
                      -0.036743_wp, -0.036935_wp, -0.037037_wp, -0.037211_wp, -0.034985_wp, &
                      -0.072742_wp, 0.161109_wp, 0.023541_wp, 0.030751_wp, 0.021879_wp, &
                      0.024292_wp, 0.021611_wp, 0.022266_wp, 0.020134_wp, 0.019873_wp, &
                      0.019154_wp, 0.019287_wp, 0.018979_wp, 0.018874_wp, 0.018837_wp, &
                      0.018891_wp, 0.018883_wp, 0.018838_wp, 0.019742_wp, 0.019769_wp, &
                      0.021921_wp, 0.021894_wp, 0.023688_wp]
         case ("resp")
            charges = [ &
                      -0.718211_wp, 0.277791_wp, -0.053003_wp, -0.001007_wp, 0.0867_wp, &
                      0.012795_wp, 0.070388_wp, -0.041717_wp, 0.005065_wp, 0.21386_wp, &
                      -0.355253_wp, 0.42339_wp, -0.017601_wp, 0.041892_wp, 0.017264_wp, &
                      0.019128_wp, -0.006342_wp, -0.016949_wp, -0.022891_wp, -0.026958_wp, &
                      -0.010557_wp, -0.006843_wp, -0.02325_wp, -0.022209_wp, -0.00327_wp, &
                      -0.004803_wp, -0.001041_wp, -0.00148_wp, -0.03402_wp, -0.034216_wp, &
                      0.075472_wp, 0.075505_wp, 0.082374_wp]
         case ("mbis")
            charges = [ &
                      -0.679703_wp, 0.102395_wp, -0.325431_wp, -0.185492_wp, -0.215333_wp, &
                      -0.215199_wp, -0.211188_wp, -0.21762_wp, -0.214782_wp, -0.16058_wp, &
                      -0.434548_wp, 0.455029_wp, 0.066055_wp, 0.092788_wp, 0.130385_wp, &
                      0.124281_wp, 0.107032_wp, 0.108045_wp, 0.109413_wp, 0.108839_wp, &
                      0.10748_wp, 0.107687_wp, 0.106857_wp, 0.106642_wp, 0.107203_wp, &
                      0.107277_wp, 0.107599_wp, 0.107485_wp, 0.101463_wp, 0.1015_wp, &
                      0.13081_wp, 0.13074_wp, 0.13289_wp]
         case ("chelpg")
            charges = [ &
                      -0.720528_wp, 0.287798_wp, -0.063495_wp, 0.002007_wp, 0.094121_wp, &
                      0.002054_wp, 0.084995_wp, -0.049461_wp, 0.001569_wp, 0.228393_wp, &
                      -0.373294_wp, 0.423855_wp, -0.019964_wp, 0.039887_wp, 0.019394_wp, &
                      0.020951_wp, -0.00669_wp, -0.017523_wp, -0.024693_wp, -0.0287_wp, &
                      -0.009181_wp, -0.005369_wp, -0.026098_wp, -0.025053_wp, -0.00191_wp, &
                      -0.003467_wp, -0.000591_wp, -0.00105_wp, -0.036638_wp, -0.036836_wp, &
                      0.079522_wp, 0.079558_wp, 0.086435_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.281165_wp, 0.039942_wp, -0.042489_wp, -0.033693_wp, -0.035761_wp, &
                      -0.03685_wp, -0.037084_wp, -0.037181_wp, -0.037349_wp, -0.035121_wp, &
                      -0.073156_wp, 0.163618_wp, 0.024714_wp, 0.030387_wp, 0.021243_wp, &
                      0.025577_wp, 0.022177_wp, 0.022375_wp, 0.020081_wp, 0.020136_wp, &
                      0.019317_wp, 0.01931_wp, 0.018981_wp, 0.018971_wp, 0.018936_wp, &
                      0.018933_wp, 0.01894_wp, 0.018935_wp, 0.019847_wp, 0.019845_wp, &
                      0.021957_wp, 0.021951_wp, 0.023693_wp]
         case ("resp")
            charges = [ &
                      -0.733753_wp, 0.283519_wp, -0.050182_wp, -0.000621_wp, 0.087574_wp, &
                      0.01447_wp, 0.072195_wp, -0.041125_wp, 0.006216_wp, 0.217103_wp, &
                      -0.358434_wp, 0.431981_wp, -0.01619_wp, 0.039846_wp, 0.014761_wp, &
                      0.020696_wp, -0.006107_wp, -0.017299_wp, -0.023566_wp, -0.027346_wp, &
                      -0.011024_wp, -0.007731_wp, -0.024213_wp, -0.022897_wp, -0.003833_wp, &
                      -0.005423_wp, -0.00162_wp, -0.001995_wp, -0.034752_wp, -0.035032_wp, &
                      0.075929_wp, 0.076006_wp, 0.082847_wp]
         case ("mbis")
            charges = [ &
                      -0.691483_wp, 0.103799_wp, -0.32561_wp, -0.185872_wp, -0.215832_wp, &
                      -0.215673_wp, -0.2117_wp, -0.21816_wp, -0.215313_wp, -0.160617_wp, &
                      -0.436629_wp, 0.462371_wp, 0.068136_wp, 0.092186_wp, 0.129149_wp, &
                      0.126654_wp, 0.108038_wp, 0.108239_wp, 0.109405_wp, 0.109387_wp, &
                      0.107877_wp, 0.10783_wp, 0.106986_wp, 0.106947_wp, 0.107513_wp, &
                      0.107481_wp, 0.10783_wp, 0.107791_wp, 0.101727_wp, 0.101712_wp, &
                      0.131285_wp, 0.131254_wp, 0.133311_wp]
         case ("chelpg")
            charges = [ &
                      -0.736062_wp, 0.293455_wp, -0.060493_wp, 0.002338_wp, 0.094938_wp, &
                      0.0039_wp, 0.086805_wp, -0.048907_wp, 0.00285_wp, 0.23162_wp, &
                      -0.376506_wp, 0.432446_wp, -0.018538_wp, 0.037846_wp, 0.016844_wp, &
                      0.022482_wp, -0.006454_wp, -0.01787_wp, -0.025365_wp, -0.029087_wp, &
                      -0.009687_wp, -0.006299_wp, -0.027076_wp, -0.025754_wp, -0.002479_wp, &
                      -0.004092_wp, -0.001202_wp, -0.001596_wp, -0.037371_wp, -0.037653_wp, &
                      0.079985_wp, 0.080064_wp, 0.086919_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (17)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070846_wp, -0.034445_wp, -0.036698_wp, -0.036037_wp, -0.035939_wp, &
                      -0.036703_wp, -0.039771_wp, 0.070114_wp, -0.168894_wp, 0.021831_wp, &
                      0.023549_wp, 0.021592_wp, 0.019941_wp, 0.019536_wp, 0.018471_wp, &
                      0.01902_wp, 0.020514_wp, 0.01912_wp, 0.017983_wp, 0.020253_wp, &
                      0.023429_wp, 0.02017_wp, 0.029276_wp, 0.028654_wp, 0.032084_wp, &
                      0.033761_wp]
         case ("resp")
            charges = [ &
                      -0.350411_wp, 0.191578_wp, 0.041801_wp, -0.063839_wp, -0.005087_wp, &
                      0.116852_wp, -0.102229_wp, 0.107601_wp, -0.275343_wp, 0.074408_wp, &
                      0.084144_wp, 0.074917_wp, -0.029582_wp, -0.030084_wp, -0.007336_wp, &
                      -0.005016_wp, 0.006059_wp, -0.001237_wp, -0.000804_wp, -0.005812_wp, &
                      0.005111_wp, -0.011168_wp, 0.018171_wp, 0.056687_wp, 0.050999_wp, &
                      0.059621_wp]
         case ("mbis")
            charges = [ &
                      -0.421283_wp, -0.160571_wp, -0.211112_wp, -0.21536_wp, -0.212038_wp, &
                      -0.189757_wp, -0.31051_wp, 0.09564_wp, -0.255947_wp, 0.127799_wp, &
                      0.12981_wp, 0.127374_wp, 0.100767_wp, 0.100053_wp, 0.105396_wp, &
                      0.10639_wp, 0.108549_wp, 0.106224_wp, 0.103802_wp, 0.108114_wp, &
                      0.120739_wp, 0.102958_wp, 0.130352_wp, 0.134835_wp, 0.083545_wp, &
                      0.084196_wp]
         case ("chelpg")
            charges = [ &
                      -0.367921_wp, 0.201447_wp, 0.047792_wp, -0.074875_wp, -0.002723_wp, &
                      0.126769_wp, -0.117392_wp, 0.11768_wp, -0.277103_wp, 0.078245_wp, &
                      0.088558_wp, 0.078748_wp, -0.031458_wp, -0.031953_wp, -0.008645_wp, &
                      -0.006317_wp, 0.007882_wp, 0.00057_wp, -0.001043_wp, -0.006097_wp, &
                      0.003462_wp, -0.012876_wp, 0.020987_wp, 0.059645_wp, 0.048787_wp, &
                      0.057829_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072962_wp, -0.035154_wp, -0.037262_wp, -0.036751_wp, -0.036003_wp, &
                      -0.037641_wp, -0.039408_wp, 0.070706_wp, -0.185876_wp, 0.022082_wp, &
                      0.023509_wp, 0.021952_wp, 0.020195_wp, 0.019999_wp, 0.019194_wp, &
                      0.01949_wp, 0.020362_wp, 0.019697_wp, 0.019823_wp, 0.020928_wp, &
                      0.019374_wp, 0.022377_wp, 0.03539_wp, 0.029284_wp, 0.037482_wp, &
                      0.039176_wp]
         case ("resp")
            charges = [ &
                      -0.366492_wp, 0.20901_wp, 0.044683_wp, -0.061137_wp, 0.006134_wp, &
                      0.115323_wp, -0.089493_wp, 0.107586_wp, -0.306601_wp, 0.07696_wp, &
                      0.086284_wp, 0.077578_wp, -0.034171_wp, -0.034156_wp, -0.008596_wp, &
                      -0.007129_wp, 0.003157_wp, -0.003369_wp, -0.002729_wp, -0.009071_wp, &
                      -0.00406_wp, -0.008807_wp, 0.026227_wp, 0.052531_wp, 0.060942_wp, &
                      0.069394_wp]
         case ("mbis")
            charges = [ &
                      -0.432236_wp, -0.160607_wp, -0.213801_wp, -0.218413_wp, -0.213753_wp, &
                      -0.191258_wp, -0.314969_wp, 0.090407_wp, -0.277478_wp, 0.130413_wp, &
                      0.131911_wp, 0.130192_wp, 0.101593_wp, 0.101273_wp, 0.1073_wp, &
                      0.107835_wp, 0.108995_wp, 0.107918_wp, 0.10734_wp, 0.109618_wp, &
                      0.113679_wp, 0.107701_wp, 0.141699_wp, 0.136149_wp, 0.093948_wp, &
                      0.094507_wp]
         case ("chelpg")
            charges = [ &
                      -0.38436_wp, 0.219223_wp, 0.050764_wp, -0.072578_wp, 0.009799_wp, &
                      0.124519_wp, -0.104206_wp, 0.117409_wp, -0.308302_wp, 0.080874_wp, &
                      0.090771_wp, 0.081486_wp, -0.036116_wp, -0.036093_wp, -0.009956_wp, &
                      -0.008485_wp, 0.005001_wp, -0.001546_wp, -0.003283_wp, -0.009667_wp, &
                      -0.005625_wp, -0.010422_wp, 0.028954_wp, 0.055393_wp, 0.058804_wp, &
                      0.06764_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073744_wp, -0.03542_wp, -0.037425_wp, -0.037037_wp, -0.035943_wp, &
                      -0.038043_wp, -0.039153_wp, 0.070961_wp, -0.192888_wp, 0.022141_wp, &
                      0.023418_wp, 0.022123_wp, 0.020186_wp, 0.020187_wp, 0.019614_wp, &
                      0.019633_wp, 0.020019_wp, 0.019976_wp, 0.020923_wp, 0.021087_wp, &
                      0.01712_wp, 0.023423_wp, 0.038374_wp, 0.029294_wp, 0.039724_wp, &
                      0.041413_wp]
         case ("resp")
            charges = [ &
                      -0.372436_wp, 0.215833_wp, 0.046199_wp, -0.060594_wp, 0.01199_wp, &
                      0.11361_wp, -0.083619_wp, 0.107552_wp, -0.319639_wp, 0.077822_wp, &
                      0.086914_wp, 0.078623_wp, -0.03618_wp, -0.035739_wp, -0.008935_wp, &
                      -0.008152_wp, 0.001559_wp, -0.003995_wp, -0.003034_wp, -0.010847_wp, &
                      -0.008887_wp, -0.007294_wp, 0.030521_wp, 0.050137_wp, 0.065134_wp, &
                      0.073456_wp]
         case ("mbis")
            charges = [ &
                      -0.436216_wp, -0.160578_wp, -0.214731_wp, -0.219469_wp, -0.214393_wp, &
                      -0.191556_wp, -0.316551_wp, 0.08832_wp, -0.286418_wp, 0.131315_wp, &
                      0.132543_wp, 0.1313_wp, 0.101709_wp, 0.101755_wp, 0.108257_wp, &
                      0.108287_wp, 0.108632_wp, 0.10867_wp, 0.109391_wp, 0.109967_wp, &
                      0.109612_wp, 0.109896_wp, 0.147134_wp, 0.13614_wp, 0.098222_wp, &
                      0.098725_wp]
         case ("chelpg")
            charges = [ &
                      -0.390463_wp, 0.226183_wp, 0.052359_wp, -0.072312_wp, 0.016359_wp, &
                      0.122401_wp, -0.098091_wp, 0.117275_wp, -0.321316_wp, 0.08177_wp, &
                      0.091434_wp, 0.082564_wp, -0.038151_wp, -0.037701_wp, -0.010323_wp, &
                      -0.009539_wp, 0.003432_wp, -0.002147_wp, -0.003752_wp, -0.011605_wp, &
                      -0.010401_wp, -0.008854_wp, 0.033197_wp, 0.052946_wp, 0.063024_wp, &
                      0.071712_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (18)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070037_wp, -0.03394_wp, -0.036174_wp, -0.035174_wp, -0.034269_wp, &
                      -0.043061_wp, 0.041528_wp, -0.250842_wp, 0.023856_wp, 0.021678_wp, &
                      0.021861_wp, 0.019342_wp, 0.019658_wp, 0.018937_wp, 0.018513_wp, &
                      0.018866_wp, 0.019949_wp, 0.019265_wp, 0.021093_wp, 0.024259_wp, &
                      0.018328_wp, 0.018107_wp, 0.031349_wp, 0.146915_wp]
         case ("resp")
            charges = [ &
                      -0.333527_wp, 0.186072_wp, -0.013783_wp, 0.001159_wp, 0.070414_wp, &
                      -0.090273_wp, 0.25095_wp, -0.632684_wp, 0.080595_wp, 0.072951_wp, &
                      0.072825_wp, -0.027621_wp, -0.028419_wp, 0.006959_wp, 0.0017_wp, &
                      -0.011404_wp, -0.006032_wp, -0.018909_wp, -0.030384_wp, 0.03331_wp, &
                      0.016703_wp, -0.026146_wp, 0.049674_wp, 0.375871_wp]
         case ("mbis")
            charges = [ &
                      -0.421722_wp, -0.159776_wp, -0.213548_wp, -0.215964_wp, -0.181894_wp, &
                      -0.322807_wp, 0.094898_wp, -0.613485_wp, 0.130549_wp, 0.127742_wp, &
                      0.128048_wp, 0.100053_wp, 0.100671_wp, 0.106903_wp, 0.106156_wp, &
                      0.106737_wp, 0.108637_wp, 0.102599_wp, 0.105597_wp, 0.134758_wp, &
                      0.113124_wp, 0.056369_wp, 0.093437_wp, 0.412926_wp]
         case ("chelpg")
            charges = [ &
                      -0.351531_wp, 0.202264_wp, -0.021359_wp, -0.001895_wp, 0.084386_wp, &
                      -0.106274_wp, 0.26239_wp, -0.635044_wp, 0.084445_wp, 0.077069_wp, &
                      0.076905_wp, -0.030528_wp, -0.031299_wp, 0.008091_wp, 0.002734_wp, &
                      -0.011007_wp, -0.005719_wp, -0.02131_wp, -0.033046_wp, 0.036411_wp, &
                      0.019451_wp, -0.02893_wp, 0.047338_wp, 0.376458_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072736_wp, -0.034822_wp, -0.036834_wp, -0.035663_wp, -0.033738_wp, &
                      -0.042621_wp, 0.040143_wp, -0.278026_wp, 0.02382_wp, 0.022024_wp, &
                      0.022064_wp, 0.020003_wp, 0.020066_wp, 0.019473_wp, 0.01938_wp, &
                      0.020079_wp, 0.020242_wp, 0.021822_wp, 0.02232_wp, 0.021664_wp, &
                      0.024659_wp, 0.024041_wp, 0.03062_wp, 0.162029_wp]
         case ("resp")
            charges = [ &
                      -0.353603_wp, 0.205429_wp, -0.006079_wp, 0.000683_wp, 0.075073_wp, &
                      -0.07533_wp, 0.284204_wp, -0.723537_wp, 0.083522_wp, 0.075911_wp, &
                      0.075792_wp, -0.031805_wp, -0.033575_wp, 0.003536_wp, -0.00031_wp, &
                      -0.013347_wp, -0.007684_wp, -0.019407_wp, -0.031608_wp, 0.021389_wp, &
                      0.023725_wp, -0.020219_wp, 0.040641_wp, 0.426599_wp]
         case ("mbis")
            charges = [ &
                      -0.435318_wp, -0.160033_wp, -0.216866_wp, -0.219713_wp, -0.184437_wp, &
                      -0.32513_wp, 0.102675_wp, -0.683548_wp, 0.133207_wp, 0.131076_wp, &
                      0.131118_wp, 0.101737_wp, 0.101886_wp, 0.108577_wp, 0.108431_wp, &
                      0.109478_wp, 0.109849_wp, 0.107216_wp, 0.107974_wp, 0.129946_wp, &
                      0.12502_wp, 0.066854_wp, 0.092556_wp, 0.457456_wp]
         case ("chelpg")
            charges = [ &
                      -0.371815_wp, 0.221556_wp, -0.012834_wp, -0.002851_wp, 0.089156_wp, &
                      -0.090801_wp, 0.295397_wp, -0.725836_wp, 0.087449_wp, 0.080063_wp, &
                      0.07991_wp, -0.034729_wp, -0.036474_wp, 0.00447_wp, 0.000533_wp, &
                      -0.012905_wp, -0.007316_wp, -0.021872_wp, -0.034316_wp, 0.024347_wp, &
                      0.026352_wp, -0.02297_wp, 0.038327_wp, 0.427161_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073028_wp, -0.034912_wp, -0.036889_wp, -0.035692_wp, -0.03364_wp, &
                      -0.042533_wp, 0.039985_wp, -0.281107_wp, 0.023783_wp, 0.022061_wp, &
                      0.022062_wp, 0.020086_wp, 0.020083_wp, 0.0195_wp, 0.019506_wp, &
                      0.020267_wp, 0.020213_wp, 0.022211_wp, 0.022391_wp, 0.021225_wp, &
                      0.025521_wp, 0.024842_wp, 0.030365_wp, 0.163708_wp]
         case ("resp")
            charges = [ &
                      -0.355684_wp, 0.207595_wp, -0.005181_wp, 0.000619_wp, 0.0757_wp, &
                      -0.073485_wp, 0.287909_wp, -0.73393_wp, 0.083758_wp, 0.076211_wp, &
                      0.076047_wp, -0.032259_wp, -0.034212_wp, 0.003081_wp, -0.000487_wp, &
                      -0.013452_wp, -0.007982_wp, -0.019278_wp, -0.031891_wp, 0.019706_wp, &
                      0.024804_wp, -0.01922_wp, 0.03928_wp, 0.432352_wp]
         case ("mbis")
            charges = [ &
                      -0.436714_wp, -0.160045_wp, -0.217185_wp, -0.220059_wp, -0.184677_wp, &
                      -0.325244_wp, 0.103619_wp, -0.691444_wp, 0.133426_wp, 0.131427_wp, &
                      0.131397_wp, 0.101936_wp, 0.101962_wp, 0.108691_wp, 0.108723_wp, &
                      0.109863_wp, 0.10985_wp, 0.107903_wp, 0.108095_wp, 0.129093_wp, &
                      0.126612_wp, 0.068276_wp, 0.092132_wp, 0.46237_wp]
         case ("chelpg")
            charges = [ &
                      -0.373913_wp, 0.223713_wp, -0.011832_wp, -0.002979_wp, 0.089789_wp, &
                      -0.08888_wp, 0.299071_wp, -0.736227_wp, 0.087692_wp, 0.080365_wp, &
                      0.080166_wp, -0.035185_wp, -0.037114_wp, 0.003989_wp, 0.000332_wp, &
                      -0.013003_wp, -0.007606_wp, -0.021748_wp, -0.034604_wp, 0.022642_wp, &
                      0.027416_wp, -0.021963_wp, 0.036965_wp, 0.432916_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (19)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.25085_wp, 0.04152_wp, -0.042987_wp, -0.033871_wp, -0.03533_wp, &
                      -0.034291_wp, -0.069873_wp, 0.146899_wp, 0.031378_wp, 0.018041_wp, &
                      0.018348_wp, 0.024278_wp, 0.019901_wp, 0.019472_wp, 0.019955_wp, &
                      0.018859_wp, 0.019327_wp, 0.021169_wp, 0.024105_wp, 0.021808_wp, &
                      0.022146_wp]
         case ("resp")
            charges = [ &
                      -0.63283_wp, 0.23668_wp, -0.012571_wp, 0.176536_wp, 0.078512_wp, &
                      -0.090396_wp, -0.362606_wp, 0.374979_wp, 0.05239_wp, -0.019787_wp, &
                      0.002204_wp, 0.018031_wp, -0.022496_wp, -0.026488_wp, -0.011647_wp, &
                      -0.014897_wp, 0.007882_wp, -0.000758_wp, 0.088963_wp, 0.079253_wp, &
                      0.079045_wp]
         case ("mbis")
            charges = [ &
                      -0.613673_wp, 0.095315_wp, -0.321973_wp, -0.162495_wp, -0.213403_wp, &
                      -0.187064_wp, -0.420825_wp, 0.413022_wp, 0.093462_wp, 0.056319_wp, &
                      0.112925_wp, 0.134645_wp, 0.101248_wp, 0.100488_wp, 0.108966_wp, &
                      0.106987_wp, 0.103119_wp, 0.106092_wp, 0.130796_wp, 0.127707_wp, &
                      0.12835_wp]
         case ("chelpg")
            charges = [ &
                      -0.634628_wp, 0.242265_wp, -0.011708_wp, 0.185151_wp, 0.084686_wp, &
                      -0.099601_wp, -0.378847_wp, 0.375156_wp, 0.051222_wp, -0.021151_wp, &
                      0.002067_wp, 0.018006_wp, -0.024207_wp, -0.028234_wp, -0.013007_wp, &
                      -0.01626_wp, 0.009637_wp, 0.000875_wp, 0.093143_wp, 0.082826_wp, &
                      0.082608_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.278327_wp, 0.040109_wp, -0.042494_wp, -0.034566_wp, -0.035755_wp, &
                      -0.03372_wp, -0.072582_wp, 0.162152_wp, 0.030631_wp, 0.023973_wp, &
                      0.02482_wp, 0.021681_wp, 0.020426_wp, 0.020347_wp, 0.020303_wp, &
                      0.020159_wp, 0.021964_wp, 0.022451_wp, 0.024009_wp, 0.022176_wp, &
                      0.02225_wp]
         case ("resp")
            charges = [ &
                      -0.725098_wp, 0.271251_wp, 0.00631_wp, 0.198929_wp, 0.083629_wp, &
                      -0.096509_wp, -0.385049_wp, 0.426162_wp, 0.043063_wp, -0.013995_wp, &
                      0.008783_wp, 0.005277_wp, -0.028002_wp, -0.030622_wp, -0.013723_wp, &
                      -0.017049_wp, 0.009406_wp, -4.9e-05_wp, 0.092069_wp, 0.083123_wp, &
                      0.082093_wp]
         case ("mbis")
            charges = [ &
                      -0.684453_wp, 0.10318_wp, -0.324221_wp, -0.162523_wp, -0.217078_wp, &
                      -0.189881_wp, -0.434541_wp, 0.457975_wp, 0.09256_wp, 0.066802_wp, &
                      0.125051_wp, 0.129791_wp, 0.102555_wp, 0.102427_wp, 0.110273_wp, &
                      0.109878_wp, 0.107884_wp, 0.108569_wp, 0.133379_wp, 0.131112_wp, &
                      0.131269_wp]
         case ("chelpg")
            charges = [ &
                      -0.72679_wp, 0.276206_wp, 0.009064_wp, 0.207851_wp, 0.090389_wp, &
                      -0.107025_wp, -0.401862_wp, 0.426285_wp, 0.041966_wp, -0.015263_wp, &
                      0.008228_wp, 0.00479_wp, -0.029771_wp, -0.032416_wp, -0.0152_wp, &
                      -0.018547_wp, 0.011336_wp, 0.001785_wp, 0.096389_wp, 0.086813_wp, &
                      0.085772_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.281118_wp, 0.039964_wp, -0.042408_wp, -0.034621_wp, -0.035775_wp, &
                      -0.033626_wp, -0.072847_wp, 0.163672_wp, 0.0304_wp, 0.024691_wp, &
                      0.025609_wp, 0.021287_wp, 0.020451_wp, 0.020461_wp, 0.020281_wp, &
                      0.020336_wp, 0.022322_wp, 0.02252_wp, 0.023967_wp, 0.022214_wp, &
                      0.022227_wp]
         case ("resp")
            charges = [ &
                      -0.734561_wp, 0.274703_wp, 0.008327_wp, 0.201164_wp, 0.084115_wp, &
                      -0.096978_wp, -0.387164_wp, 0.431369_wp, 0.041813_wp, -0.013111_wp, &
                      0.009718_wp, 0.003692_wp, -0.02861_wp, -0.030982_wp, -0.014034_wp, &
                      -0.017166_wp, 0.009721_wp, -0.000115_wp, 0.092295_wp, 0.083492_wp, &
                      0.082313_wp]
         case ("mbis")
            charges = [ &
                      -0.691605_wp, 0.104035_wp, -0.324316_wp, -0.162508_wp, -0.217377_wp, &
                      -0.190119_wp, -0.43581_wp, 0.462424_wp, 0.092175_wp, 0.068078_wp, &
                      0.126505_wp, 0.129021_wp, 0.102627_wp, 0.102661_wp, 0.110282_wp, &
                      0.110237_wp, 0.108517_wp, 0.108688_wp, 0.133564_wp, 0.131441_wp, &
                      0.131486_wp]
         case ("chelpg")
            charges = [ &
                      -0.736248_wp, 0.279595_wp, 0.011278_wp, 0.210122_wp, 0.090926_wp, &
                      -0.107629_wp, -0.404028_wp, 0.431491_wp, 0.040719_wp, -0.014366_wp, &
                      0.009122_wp, 0.003155_wp, -0.030387_wp, -0.032783_wp, -0.015522_wp, &
                      -0.018675_wp, 0.011669_wp, 0.00174_wp, 0.096627_wp, 0.087192_wp, &
                      0.086001_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (20)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.035761_wp, 0.025553_wp, 0.025439_wp, -0.024126_wp, 0.028894_wp, &
                      -0.0871_wp, 0.029456_wp, 0.027783_wp, -0.034859_wp, 0.019842_wp, &
                      0.019838_wp, -0.033659_wp, 0.020261_wp, 0.019634_wp, -0.06972_wp, &
                      0.024015_wp, 0.022058_wp, 0.022457_wp]
         case ("resp")
            charges = [ &
                      -0.05145_wp, 0.049028_wp, 0.033302_wp, -0.053382_wp, 0.109764_wp, &
                      -0.529241_wp, 0.19075_wp, 0.206596_wp, -0.005191_wp, 0.008635_wp, &
                      0.017796_wp, 0.191353_wp, -0.033355_wp, -0.029664_wp, -0.320822_wp, &
                      0.077551_wp, 0.068612_wp, 0.069717_wp]
         case ("mbis")
            charges = [ &
                      -0.208496_wp, 0.119756_wp, 0.114995_wp, -0.065704_wp, 0.120945_wp, &
                      -0.433495_wp, 0.170506_wp, 0.166879_wp, -0.214945_wp, 0.108543_wp, &
                      0.115052_wp, -0.163461_wp, 0.102013_wp, 0.10041_wp, -0.420354_wp, &
                      0.13051_wp, 0.128088_wp, 0.128763_wp]
         case ("chelpg")
            charges = [ &
                      -0.057051_wp, 0.050135_wp, 0.034287_wp, -0.04771_wp, 0.10887_wp, &
                      -0.534159_wp, 0.191904_wp, 0.208095_wp, -0.00795_wp, 0.009018_wp, &
                      0.01836_wp, 0.204428_wp, -0.035748_wp, -0.031997_wp, -0.33796_wp, &
                      0.081496_wp, 0.072453_wp, 0.073528_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036599_wp, 0.02678_wp, 0.027043_wp, -0.026382_wp, 0.031155_wp, &
                      -0.09321_wp, 0.030379_wp, 0.02932_wp, -0.035015_wp, 0.021087_wp, &
                      0.019354_wp, -0.033721_wp, 0.020944_wp, 0.020522_wp, -0.071018_wp, &
                      0.024282_wp, 0.022407_wp, 0.022673_wp]
         case ("resp")
            charges = [ &
                      -0.050997_wp, 0.050608_wp, 0.035356_wp, -0.062403_wp, 0.11778_wp, &
                      -0.554982_wp, 0.198899_wp, 0.216905_wp, 0.001067_wp, 0.008697_wp, &
                      0.014486_wp, 0.201584_wp, -0.035302_wp, -0.031103_wp, -0.332965_wp, &
                      0.07985_wp, 0.070676_wp, 0.071844_wp]
         case ("mbis")
            charges = [ &
                      -0.213046_wp, 0.123114_wp, 0.118767_wp, -0.069436_wp, 0.127129_wp, &
                      -0.455236_wp, 0.177163_wp, 0.17455_wp, -0.215541_wp, 0.111103_wp, &
                      0.1144_wp, -0.163418_wp, 0.1033_wp, 0.102058_wp, -0.427964_wp, &
                      0.132404_wp, 0.130096_wp, 0.130563_wp]
         case ("chelpg")
            charges = [ &
                      -0.056712_wp, 0.051726_wp, 0.036366_wp, -0.056868_wp, 0.116934_wp, &
                      -0.559762_wp, 0.200014_wp, 0.21837_wp, -0.001268_wp, 0.00897_wp, &
                      0.014936_wp, 0.214585_wp, -0.037705_wp, -0.033446_wp, -0.350147_wp, &
                      0.083821_wp, 0.074524_wp, 0.075661_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037371_wp, 0.027993_wp, 0.028996_wp, -0.028511_wp, 0.033629_wp, &
                      -0.099899_wp, 0.031028_wp, 0.030731_wp, -0.035085_wp, 0.022569_wp, &
                      0.018524_wp, -0.033644_wp, 0.021667_wp, 0.021691_wp, -0.072296_wp, &
                      0.024492_wp, 0.022767_wp, 0.022724_wp]
         case ("resp")
            charges = [ &
                      -0.049996_wp, 0.051956_wp, 0.037957_wp, -0.071872_wp, 0.126403_wp, &
                      -0.582181_wp, 0.206994_wp, 0.227424_wp, 0.007794_wp, 0.009052_wp, &
                      0.010179_wp, 0.21318_wp, -0.037562_wp, -0.032366_wp, -0.345617_wp, &
                      0.082091_wp, 0.072821_wp, 0.073745_wp]
         case ("mbis")
            charges = [ &
                      -0.217887_wp, 0.126465_wp, 0.123167_wp, -0.072657_wp, 0.133665_wp, &
                      -0.478521_wp, 0.183716_wp, 0.18236_wp, -0.215724_wp, 0.114028_wp, &
                      0.113012_wp, -0.163272_wp, 0.104584_wp, 0.104146_wp, -0.435394_wp, &
                      0.134173_wp, 0.132101_wp, 0.132042_wp]
         case ("chelpg")
            charges = [ &
                      -0.055809_wp, 0.053076_wp, 0.038978_wp, -0.06649_wp, 0.125595_wp, &
                      -0.586788_wp, 0.208062_wp, 0.228851_wp, 0.005943_wp, 0.009195_wp, &
                      0.010497_wp, 0.226109_wp, -0.039978_wp, -0.034727_wp, -0.362833_wp, &
                      0.086085_wp, 0.076669_wp, 0.077564_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (21)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.01586_wp, -0.047656_wp, -0.029145_wp, 0.039828_wp, 0.039839_wp, &
                      -0.131214_wp, 0.025199_wp, 0.079853_wp, 0.023689_wp, -0.032896_wp, &
                      0.025404_wp, 0.025415_wp, -0.068429_wp, 0.021146_wp, 0.021138_wp, &
                      0.02369_wp]
         case ("resp")
            charges = [ &
                      -0.140899_wp, 0.124321_wp, 0.03778_wp, 0.078452_wp, 0.079601_wp, &
                      -0.57787_wp, 0.086196_wp, 0.338194_wp, 0.074331_wp, 0.162854_wp, &
                      0.020754_wp, 0.019141_wp, -0.33555_wp, -0.021884_wp, -0.02105_wp, &
                      0.07563_wp]
         case ("mbis")
            charges = [ &
                      -0.198312_wp, 0.053626_wp, -0.215337_wp, 0.136163_wp, 0.136174_wp, &
                      -0.442759_wp, 0.131815_wp, 0.274666_wp, 0.130046_wp, -0.168325_wp, &
                      0.122514_wp, 0.122438_wp, -0.418589_wp, 0.102918_wp, 0.102899_wp, &
                      0.130065_wp]
         case ("chelpg")
            charges = [ &
                      -0.172467_wp, 0.143684_wp, 0.048491_wp, 0.084657_wp, 0.085804_wp, &
                      -0.589302_wp, 0.091305_wp, 0.341438_wp, 0.078852_wp, 0.174344_wp, &
                      0.018829_wp, 0.017189_wp, -0.356417_wp, -0.023684_wp, -0.022872_wp, &
                      0.080146_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.016497_wp, -0.053909_wp, -0.029146_wp, 0.044157_wp, 0.044169_wp, &
                      -0.143955_wp, 0.025572_wp, 0.086577_wp, 0.023974_wp, -0.032136_wp, &
                      0.025539_wp, 0.025562_wp, -0.069826_wp, 0.022977_wp, 0.022972_wp, &
                      0.023973_wp]
         case ("resp")
            charges = [ &
                      -0.152249_wp, 0.116271_wp, 0.049138_wp, 0.087823_wp, 0.088974_wp, &
                      -0.606616_wp, 0.09023_wp, 0.357191_wp, 0.078446_wp, 0.180124_wp, &
                      0.017518_wp, 0.015871_wp, -0.357336_wp, -0.023001_wp, -0.022144_wp, &
                      0.07976_wp]
         case ("mbis")
            charges = [ &
                      -0.209492_wp, 0.051852_wp, -0.214414_wp, 0.145786_wp, 0.145796_wp, &
                      -0.473661_wp, 0.134137_wp, 0.293235_wp, 0.132195_wp, -0.166998_wp, &
                      0.122912_wp, 0.122856_wp, -0.427939_wp, 0.105769_wp, 0.105758_wp, &
                      0.13221_wp]
         case ("chelpg")
            charges = [ &
                      -0.184964_wp, 0.136132_wp, 0.060879_wp, 0.094245_wp, 0.095394_wp, &
                      -0.618275_wp, 0.095433_wp, 0.360487_wp, 0.083017_wp, 0.191574_wp, &
                      0.015357_wp, 0.013682_wp, -0.378485_wp, -0.024819_wp, -0.023984_wp, &
                      0.084327_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.016863_wp, -0.058307_wp, -0.02911_wp, 0.04746_wp, 0.047472_wp, &
                      -0.153162_wp, 0.025824_wp, 0.090971_wp, 0.024003_wp, -0.031362_wp, &
                      0.025327_wp, 0.02536_wp, -0.070756_wp, 0.024573_wp, 0.024571_wp, &
                      0.024_wp]
         case ("resp")
            charges = [ &
                      -0.159004_wp, 0.109062_wp, 0.055961_wp, 0.094776_wp, 0.095927_wp, &
                      -0.625912_wp, 0.093124_wp, 0.369755_wp, 0.08121_wp, 0.194484_wp, &
                      0.014805_wp, 0.013138_wp, -0.373426_wp, -0.023655_wp, -0.022781_wp, &
                      0.082537_wp]
         case ("mbis")
            charges = [ &
                      -0.217579_wp, 0.050422_wp, -0.21312_wp, 0.152996_wp, 0.153005_wp, &
                      -0.495313_wp, 0.135684_wp, 0.305688_wp, 0.133329_wp, -0.16596_wp, &
                      0.122565_wp, 0.122527_wp, -0.434046_wp, 0.108234_wp, 0.10823_wp, &
                      0.133341_wp]
         case ("chelpg")
            charges = [ &
                      -0.192321_wp, 0.129169_wp, 0.068254_wp, 0.101313_wp, 0.102463_wp, &
                      -0.637674_wp, 0.098382_wp, 0.373073_wp, 0.085812_wp, 0.205927_wp, &
                      0.012515_wp, 0.010819_wp, -0.394744_wp, -0.025487_wp, -0.024636_wp, &
                      0.087135_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (22)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06751_wp, 0.025079_wp, 0.026594_wp, 0.025069_wp, -0.031419_wp, &
                      0.022462_wp, 0.022556_wp, -0.034871_wp, 0.028331_wp, 0.028304_wp, &
                      0.002256_wp, -0.125479_wp, 0.03927_wp, 0.039354_wp]
         case ("resp")
            charges = [ &
                      -0.368119_wp, 0.081756_wp, 0.096257_wp, 0.081119_wp, 0.183269_wp, &
                      -0.022797_wp, -0.022656_wp, 0.147404_wp, 0.010431_wp, 0.008387_wp, &
                      -0.385307_wp, -0.136295_wp, 0.162823_wp, 0.163728_wp]
         case ("mbis")
            charges = [ &
                      -0.425017_wp, 0.132825_wp, 0.134765_wp, 0.132675_wp, -0.162124_wp, &
                      0.10599_wp, 0.105892_wp, -0.202655_wp, 0.133756_wp, 0.133821_wp, &
                      -0.290955_wp, -0.138852_wp, 0.169916_wp, 0.169959_wp]
         case ("chelpg")
            charges = [ &
                      -0.383401_wp, 0.085074_wp, 0.100198_wp, 0.084428_wp, 0.191632_wp, &
                      -0.024577_wp, -0.024419_wp, 0.153653_wp, 0.009143_wp, 0.007079_wp, &
                      -0.395322_wp, -0.134942_wp, 0.16528_wp, 0.166171_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06913_wp, 0.025408_wp, 0.026987_wp, 0.025412_wp, -0.029423_wp, &
                      0.025609_wp, 0.025686_wp, -0.033889_wp, 0.029368_wp, 0.029314_wp, &
                      0.010176_wp, -0.159697_wp, 0.047039_wp, 0.047136_wp]
         case ("resp")
            charges = [ &
                      -0.386625_wp, 0.084791_wp, 0.098694_wp, 0.084209_wp, 0.205555_wp, &
                      -0.023848_wp, -0.023706_wp, 0.139759_wp, 0.011431_wp, 0.009392_wp, &
                      -0.347951_wp, -0.189694_wp, 0.168549_wp, 0.169442_wp]
         case ("mbis")
            charges = [ &
                      -0.437084_wp, 0.135519_wp, 0.137592_wp, 0.135394_wp, -0.159341_wp, &
                      0.110507_wp, 0.110386_wp, -0.207468_wp, 0.135844_wp, 0.135863_wp, &
                      -0.275907_wp, -0.182393_wp, 0.180506_wp, 0.180577_wp]
         case ("chelpg")
            charges = [ &
                      -0.402114_wp, 0.088155_wp, 0.102672_wp, 0.087564_wp, 0.214201_wp, &
                      -0.025683_wp, -0.025523_wp, 0.145732_wp, 0.010197_wp, 0.008139_wp, &
                      -0.357808_wp, -0.188355_wp, 0.170972_wp, 0.171851_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069463_wp, 0.025382_wp, 0.027048_wp, 0.025391_wp, -0.028807_wp, &
                      0.026503_wp, 0.026573_wp, -0.033592_wp, 0.029506_wp, 0.029446_wp, &
                      0.012412_wp, -0.168615_wp, 0.049057_wp, 0.049157_wp]
         case ("resp")
            charges = [ &
                      -0.389841_wp, 0.085113_wp, 0.098905_wp, 0.08455_wp, 0.211302_wp, &
                      -0.023973_wp, -0.023835_wp, 0.134707_wp, 0.01204_wp, 0.01001_wp, &
                      -0.333506_wp, -0.204608_wp, 0.169123_wp, 0.170012_wp]
         case ("mbis")
            charges = [ &
                      -0.439529_wp, 0.135908_wp, 0.138143_wp, 0.135792_wp, -0.158618_wp, &
                      0.111784_wp, 0.111652_wp, -0.208521_wp, 0.136111_wp, 0.136119_wp, &
                      -0.271314_wp, -0.193812_wp, 0.183102_wp, 0.18318_wp]
         case ("chelpg")
            charges = [ &
                      -0.405371_wp, 0.088487_wp, 0.102888_wp, 0.087914_wp, 0.22003_wp, &
                      -0.025821_wp, -0.025666_wp, 0.140543_wp, 0.010836_wp, 0.008786_wp, &
                      -0.343289_wp, -0.203276_wp, 0.171532_wp, 0.172406_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (23)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071133_wp, -0.035555_wp, -0.036744_wp, -0.037463_wp, -0.037376_wp, &
                      -0.036765_wp, -0.0374_wp, -0.037615_wp, -0.036617_wp, -0.036534_wp, &
                      -0.037041_wp, -0.036959_wp, -0.03507_wp, -0.034172_wp, -0.036346_wp, &
                      0.00084_wp, -0.124738_wp, 0.021492_wp, 0.023172_wp, 0.021835_wp, &
                      0.019237_wp, 0.019847_wp, 0.018011_wp, 0.019426_wp, 0.019296_wp, &
                      0.01873_wp, 0.018741_wp, 0.019352_wp, 0.019138_wp, 0.01781_wp, &
                      0.019093_wp, 0.018691_wp, 0.01853_wp, 0.018427_wp, 0.02003_wp, &
                      0.018468_wp, 0.018478_wp, 0.020023_wp, 0.019297_wp, 0.019298_wp, &
                      0.019482_wp, 0.020514_wp, 0.020725_wp, 0.021785_wp, 0.022294_wp, &
                      0.021309_wp, 0.027803_wp, 0.027921_wp, 0.039027_wp, 0.039399_wp]
         case ("resp")
            charges = [ &
                      -0.420856_wp, 0.289803_wp, -0.142677_wp, -0.087536_wp, 0.045036_wp, &
                      -0.135941_wp, 0.013661_wp, 0.019779_wp, -0.06364_wp, -0.042444_wp, &
                      -0.016266_wp, -0.001856_wp, -0.106482_wp, -0.018496_wp, 0.21207_wp, &
                      -0.415394_wp, -0.145402_wp, 0.095367_wp, 0.092594_wp, 0.092144_wp, &
                      -0.037472_wp, -0.042589_wp, 0.037274_wp, 0.034385_wp, 0.044824_wp, &
                      0.026919_wp, 0.003815_wp, 0.010527_wp, 0.028079_wp, 0.029879_wp, &
                      0.008386_wp, 0.007968_wp, 0.006109_wp, 0.013363_wp, 0.01656_wp, &
                      0.016516_wp, 0.013553_wp, 0.019704_wp, 0.025752_wp, 0.009525_wp, &
                      0.021092_wp, -6.3e-05_wp, 0.023793_wp, 0.035569_wp, 0.007251_wp, &
                      0.013624_wp, 0.001403_wp, 0.010059_wp, 0.177784_wp, 0.17295_wp]
         case ("mbis")
            charges = [ &
                      -0.426318_wp, -0.152295_wp, -0.219349_wp, -0.208158_wp, -0.203992_wp, &
                      -0.217099_wp, -0.206398_wp, -0.205658_wp, -0.212136_wp, -0.211883_wp, &
                      -0.205783_wp, -0.208391_wp, -0.216962_wp, -0.212405_wp, -0.195094_wp, &
                      -0.292071_wp, -0.139554_wp, 0.127507_wp, 0.129808_wp, 0.128856_wp, &
                      0.097936_wp, 0.102114_wp, 0.1038_wp, 0.107014_wp, 0.104389_wp, &
                      0.10481_wp, 0.104155_wp, 0.103796_wp, 0.10528_wp, 0.10316_wp, &
                      0.107725_wp, 0.104384_wp, 0.105433_wp, 0.103938_wp, 0.105926_wp, &
                      0.103853_wp, 0.103798_wp, 0.105853_wp, 0.10539_wp, 0.106766_wp, &
                      0.105805_wp, 0.108384_wp, 0.107741_wp, 0.109868_wp, 0.110868_wp, &
                      0.108935_wp, 0.134947_wp, 0.131377_wp, 0.169261_wp, 0.170661_wp]
         case ("chelpg")
            charges = [ &
                      -0.438292_wp, 0.306682_wp, -0.156706_wp, -0.093836_wp, 0.060162_wp, &
                      -0.154916_wp, 0.017932_wp, 0.026991_wp, -0.073016_wp, -0.041204_wp, &
                      -0.018121_wp, 0.000624_wp, -0.113367_wp, -0.021894_wp, 0.224013_wp, &
                      -0.428355_wp, -0.143998_wp, 0.099169_wp, 0.096424_wp, 0.096158_wp, &
                      -0.040207_wp, -0.045841_wp, 0.040441_wp, 0.038159_wp, 0.047393_wp, &
                      0.02916_wp, 0.001528_wp, 0.008345_wp, 0.032081_wp, 0.033599_wp, &
                      0.007177_wp, 0.007455_wp, 0.00492_wp, 0.012547_wp, 0.01868_wp, &
                      0.018441_wp, 0.01358_wp, 0.020377_wp, 0.026609_wp, 0.009812_wp, &
                      0.021149_wp, -0.000764_wp, 0.025436_wp, 0.038041_wp, 0.008103_wp, &
                      0.014419_wp, -0.000504_wp, 0.008148_wp, 0.181051_wp, 0.176214_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073214_wp, -0.036049_wp, -0.037548_wp, -0.038061_wp, -0.03805_wp, &
                      -0.037685_wp, -0.038117_wp, -0.038248_wp, -0.037602_wp, -0.037477_wp, &
                      -0.037697_wp, -0.037358_wp, -0.035213_wp, -0.03281_wp, -0.035528_wp, &
                      0.006981_wp, -0.1523_wp, 0.021771_wp, 0.023295_wp, 0.022046_wp, &
                      0.019854_wp, 0.020489_wp, 0.018436_wp, 0.019602_wp, 0.019706_wp, &
                      0.019351_wp, 0.019294_wp, 0.019637_wp, 0.019304_wp, 0.018107_wp, &
                      0.019461_wp, 0.018971_wp, 0.019258_wp, 0.018848_wp, 0.019826_wp, &
                      0.018613_wp, 0.01869_wp, 0.019897_wp, 0.019568_wp, 0.019873_wp, &
                      0.020156_wp, 0.021137_wp, 0.021413_wp, 0.022331_wp, 0.024797_wp, &
                      0.023735_wp, 0.028505_wp, 0.028994_wp, 0.045372_wp, 0.045635_wp]
         case ("resp")
            charges = [ &
                      -0.441445_wp, 0.311397_wp, -0.152504_wp, -0.090082_wp, 0.049294_wp, &
                      -0.143989_wp, 0.017488_wp, 0.022753_wp, -0.065068_wp, -0.040687_wp, &
                      -0.017648_wp, 0.001932_wp, -0.107712_wp, -0.021651_wp, 0.214092_wp, &
                      -0.388522_wp, -0.188384_wp, 0.09928_wp, 0.095737_wp, 0.095708_wp, &
                      -0.040812_wp, -0.046355_wp, 0.039251_wp, 0.035612_wp, 0.046795_wp, &
                      0.02804_wp, 0.00379_wp, 0.010393_wp, 0.028794_wp, 0.030853_wp, &
                      0.007719_wp, 0.00701_wp, 0.006006_wp, 0.013552_wp, 0.01502_wp, &
                      0.016086_wp, 0.012434_wp, 0.017972_wp, 0.026119_wp, 0.009803_wp, &
                      0.021283_wp, -0.000384_wp, 0.024533_wp, 0.03546_wp, 0.010529_wp, &
                      0.017454_wp, 0.000452_wp, 0.011025_wp, 0.183342_wp, 0.178235_wp]
         case ("mbis")
            charges = [ &
                      -0.437407_wp, -0.151518_wp, -0.223064_wp, -0.21052_wp, -0.206168_wp, &
                      -0.220819_wp, -0.208462_wp, -0.208059_wp, -0.214957_wp, -0.214828_wp, &
                      -0.20821_wp, -0.210196_wp, -0.220194_wp, -0.213279_wp, -0.198142_wp, &
                      -0.280647_wp, -0.174722_wp, 0.130169_wp, 0.132184_wp, 0.131379_wp, &
                      0.099233_wp, 0.103466_wp, 0.105409_wp, 0.108135_wp, 0.105717_wp, &
                      0.10652_wp, 0.10573_wp, 0.104867_wp, 0.106422_wp, 0.10457_wp, &
                      0.108972_wp, 0.105442_wp, 0.107389_wp, 0.105287_wp, 0.106315_wp, &
                      0.104922_wp, 0.104995_wp, 0.106373_wp, 0.106442_wp, 0.108414_wp, &
                      0.107359_wp, 0.109886_wp, 0.109493_wp, 0.111307_wp, 0.114826_wp, &
                      0.112888_wp, 0.136265_wp, 0.133386_wp, 0.178123_wp, 0.1793_wp]
         case ("chelpg")
            charges = [ &
                      -0.459028_wp, 0.328528_wp, -0.166875_wp, -0.096613_wp, 0.065131_wp, &
                      -0.163804_wp, 0.022212_wp, 0.030282_wp, -0.07487_wp, -0.03881_wp, &
                      -0.020153_wp, 0.004892_wp, -0.114515_wp, -0.025646_wp, 0.226213_wp, &
                      -0.401465_wp, -0.186994_wp, 0.103107_wp, 0.099596_wp, 0.099756_wp, &
                      -0.043588_wp, -0.049669_wp, 0.042496_wp, 0.039493_wp, 0.049426_wp, &
                      0.030353_wp, 0.00138_wp, 0.008096_wp, 0.032947_wp, 0.03472_wp, &
                      0.006392_wp, 0.006393_wp, 0.004742_wp, 0.012658_wp, 0.017177_wp, &
                      0.018076_wp, 0.012332_wp, 0.018531_wp, 0.027076_wp, 0.010196_wp, &
                      0.021246_wp, -0.001181_wp, 0.026179_wp, 0.037934_wp, 0.011546_wp, &
                      0.018392_wp, -0.001486_wp, 0.009102_wp, 0.1866_wp, 0.181499_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074095_wp, -0.036228_wp, -0.037878_wp, -0.038294_wp, -0.038318_wp, &
                      -0.038058_wp, -0.03842_wp, -0.038492_wp, -0.038053_wp, -0.037904_wp, &
                      -0.037971_wp, -0.037505_wp, -0.035189_wp, -0.031916_wp, -0.035032_wp, &
                      0.010531_wp, -0.166951_wp, 0.02185_wp, 0.023306_wp, 0.0221_wp, &
                      0.020142_wp, 0.020793_wp, 0.018607_wp, 0.019635_wp, 0.019876_wp, &
                      0.019618_wp, 0.019546_wp, 0.019748_wp, 0.019378_wp, 0.018228_wp, &
                      0.019572_wp, 0.019065_wp, 0.019631_wp, 0.019049_wp, 0.019601_wp, &
                      0.018614_wp, 0.018768_wp, 0.019706_wp, 0.019617_wp, 0.020162_wp, &
                      0.020524_wp, 0.021284_wp, 0.021654_wp, 0.022614_wp, 0.026272_wp, &
                      0.025134_wp, 0.028668_wp, 0.029368_wp, 0.048717_wp, 0.048919_wp]
         case ("resp")
            charges = [ &
                      -0.450139_wp, 0.321088_wp, -0.156702_wp, -0.090947_wp, 0.05121_wp, &
                      -0.147258_wp, 0.019088_wp, 0.024323_wp, -0.065591_wp, -0.039258_wp, &
                      -0.018602_wp, 0.00354_wp, -0.106259_wp, -0.02338_wp, 0.210043_wp, &
                      -0.365555_wp, -0.212993_wp, 0.100833_wp, 0.09696_wp, 0.09712_wp, &
                      -0.042332_wp, -0.048034_wp, 0.040037_wp, 0.036008_wp, 0.047565_wp, &
                      0.028461_wp, 0.003787_wp, 0.010254_wp, 0.029048_wp, 0.031191_wp, &
                      0.007344_wp, 0.006513_wp, 0.006011_wp, 0.013602_wp, 0.013992_wp, &
                      0.015773_wp, 0.011765_wp, 0.016756_wp, 0.026084_wp, 0.010042_wp, &
                      0.02144_wp, -0.00062_wp, 0.024405_wp, 0.034955_wp, 0.012556_wp, &
                      0.01972_wp, 0.000447_wp, 0.012037_wp, 0.184459_wp, 0.179212_wp]
         case ("mbis")
            charges = [ &
                      -0.442006_wp, -0.151148_wp, -0.224576_wp, -0.211458_wp, -0.207061_wp, &
                      -0.222345_wp, -0.209278_wp, -0.209062_wp, -0.216063_wp, -0.216026_wp, &
                      -0.209212_wp, -0.210855_wp, -0.221491_wp, -0.213466_wp, -0.199545_wp, &
                      -0.273415_wp, -0.193535_wp, 0.131221_wp, 0.133109_wp, 0.132378_wp, &
                      0.099815_wp, 0.104081_wp, 0.106064_wp, 0.108519_wp, 0.106258_wp, &
                      0.107239_wp, 0.106419_wp, 0.10529_wp, 0.106899_wp, 0.105145_wp, &
                      0.109409_wp, 0.105837_wp, 0.108329_wp, 0.105892_wp, 0.106226_wp, &
                      0.105274_wp, 0.105493_wp, 0.106347_wp, 0.106762_wp, 0.109194_wp, &
                      0.108144_wp, 0.110276_wp, 0.110099_wp, 0.11196_wp, 0.117091_wp, &
                      0.115082_wp, 0.13653_wp, 0.134046_wp, 0.182525_wp, 0.18358_wp]
         case ("chelpg")
            charges = [ &
                      -0.467778_wp, 0.338321_wp, -0.171208_wp, -0.097573_wp, 0.067354_wp, &
                      -0.167437_wp, 0.024001_wp, 0.032024_wp, -0.075629_wp, -0.036984_wp, &
                      -0.021472_wp, 0.006688_wp, -0.112895_wp, -0.027744_wp, 0.222219_wp, &
                      -0.378452_wp, -0.211615_wp, 0.10467_wp, 0.10083_wp, 0.101182_wp, &
                      -0.045125_wp, -0.051373_wp, 0.043313_wp, 0.039932_wp, 0.050221_wp, &
                      0.030803_wp, 0.001324_wp, 0.007906_wp, 0.033265_wp, 0.035121_wp, &
                      0.005968_wp, 0.005852_wp, 0.004705_wp, 0.012663_wp, 0.016167_wp, &
                      0.017797_wp, 0.011582_wp, 0.01724_wp, 0.027095_wp, 0.010496_wp, &
                      0.021363_wp, -0.001454_wp, 0.026028_wp, 0.037401_wp, 0.01367_wp, &
                      0.020743_wp, -0.001496_wp, 0.010117_wp, 0.187706_wp, 0.18247_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (24)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070007_wp, -0.032626_wp, -0.034397_wp, -0.035381_wp, 0.00061_wp, &
                      -0.123995_wp, 0.023099_wp, 0.025594_wp, 0.023115_wp, 0.022497_wp, &
                      0.02251_wp, 0.021761_wp, 0.021773_wp, 0.028378_wp, 0.028426_wp, &
                      0.039339_wp, 0.039304_wp]
         case ("resp")
            charges = [ &
                      -0.294889_wp, 0.179097_wp, -0.017447_wp, 0.115136_wp, -0.288712_wp, &
                      -0.149031_wp, 0.063982_wp, 0.074519_wp, 0.064574_wp, -0.030146_wp, &
                      -0.030265_wp, 0.00807_wp, 0.00901_wp, 0.005679_wp, 0.007475_wp, &
                      0.141389_wp, 0.141558_wp]
         case ("mbis")
            charges = [ &
                      -0.420595_wp, -0.165272_wp, -0.211653_wp, -0.207326_wp, -0.286332_wp, &
                      -0.138625_wp, 0.12912_wp, 0.132621_wp, 0.129128_wp, 0.104794_wp, &
                      0.104825_wp, 0.111439_wp, 0.111405_wp, 0.13401_wp, 0.133932_wp, &
                      0.169279_wp, 0.169249_wp]
         case ("chelpg")
            charges = [ &
                      -0.308655_wp, 0.193287_wp, -0.028844_wp, 0.127043_wp, -0.298267_wp, &
                      -0.148188_wp, 0.06718_wp, 0.077405_wp, 0.06777_wp, -0.033002_wp, &
                      -0.033118_wp, 0.009889_wp, 0.010829_wp, 0.003408_wp, 0.005214_wp, &
                      0.143936_wp, 0.144113_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072133_wp, -0.032632_wp, -0.03246_wp, -0.034433_wp, 0.008121_wp, &
                      -0.157008_wp, 0.023457_wp, 0.025092_wp, 0.023482_wp, 0.023224_wp, &
                      0.023216_wp, 0.024864_wp, 0.024874_wp, 0.029311_wp, 0.029372_wp, &
                      0.04684_wp, 0.046815_wp]
         case ("resp")
            charges = [ &
                      -0.316361_wp, 0.199856_wp, -0.012289_wp, 0.106825_wp, -0.24419_wp, &
                      -0.202038_wp, 0.06756_wp, 0.076569_wp, 0.068159_wp, -0.034239_wp, &
                      -0.034355_wp, 0.010005_wp, 0.010952_wp, 0.006077_wp, 0.007893_wp, &
                      0.144692_wp, 0.144884_wp]
         case ("mbis")
            charges = [ &
                      -0.432499_wp, -0.165189_wp, -0.21175_wp, -0.21182_wp, -0.271905_wp, &
                      -0.180803_wp, 0.132076_wp, 0.134119_wp, 0.132102_wp, 0.106193_wp, &
                      0.106184_wp, 0.116248_wp, 0.11621_wp, 0.135881_wp, 0.135821_wp, &
                      0.179571_wp, 0.179561_wp]
         case ("chelpg")
            charges = [ &
                      -0.330441_wp, 0.214226_wp, -0.023298_wp, 0.118209_wp, -0.253416_wp, &
                      -0.201225_wp, 0.070819_wp, 0.07954_wp, 0.071416_wp, -0.037146_wp, &
                      -0.037259_wp, 0.011733_wp, 0.01268_wp, 0.003903_wp, 0.005728_wp, &
                      0.147166_wp, 0.147366_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072617_wp, -0.032608_wp, -0.031809_wp, -0.034116_wp, 0.010455_wp, &
                      -0.16647_wp, 0.023557_wp, 0.024806_wp, 0.023584_wp, 0.023327_wp, &
                      0.023315_wp, 0.025833_wp, 0.025843_wp, 0.02945_wp, 0.029514_wp, &
                      0.048981_wp, 0.048957_wp]
         case ("resp")
            charges = [ &
                      -0.322098_wp, 0.205807_wp, -0.010578_wp, 0.101349_wp, -0.226667_wp, &
                      -0.218272_wp, 0.068586_wp, 0.0768_wp, 0.069182_wp, -0.035502_wp, &
                      -0.03561_wp, 0.010717_wp, 0.011666_wp, 0.006596_wp, 0.008414_wp, &
                      0.144707_wp, 0.144902_wp]
         case ("mbis")
            charges = [ &
                      -0.435159_wp, -0.165074_wp, -0.211633_wp, -0.212907_wp, -0.267069_wp, &
                      -0.192965_wp, 0.132788_wp, 0.134154_wp, 0.132817_wp, 0.106381_wp, &
                      0.106362_wp, 0.117717_wp, 0.117678_wp, 0.136141_wp, 0.136087_wp, &
                      0.182345_wp, 0.182338_wp]
         case ("chelpg")
            charges = [ &
                      -0.336273_wp, 0.220202_wp, -0.021393_wp, 0.112433_wp, -0.235717_wp, &
                      -0.217474_wp, 0.071861_wp, 0.079799_wp, 0.072456_wp, -0.038418_wp, &
                      -0.038522_wp, 0.012407_wp, 0.013357_wp, 0.004483_wp, 0.00631_wp, &
                      0.147143_wp, 0.147346_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (25)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.065962_wp, 0.025915_wp, 0.030385_wp, 0.025918_wp, -0.031828_wp, &
                      0.029694_wp, 0.029689_wp, 0.002114_wp, -0.124848_wp, 0.039462_wp, &
                      0.039463_wp]
         case ("resp")
            charges = [ &
                      -0.338607_wp, 0.084009_wp, 0.078516_wp, 0.083876_wp, 0.322013_wp, &
                      -0.02079_wp, -0.019292_wp, -0.369217_wp, -0.149601_wp, 0.164731_wp, &
                      0.164361_wp]
         case ("mbis")
            charges = [ &
                      -0.421187_wp, 0.134148_wp, 0.139072_wp, 0.134146_wp, -0.156157_wp, &
                      0.12879_wp, 0.128805_wp, -0.285653_wp, -0.14097_wp, 0.169502_wp, &
                      0.169506_wp]
         case ("chelpg")
            charges = [ &
                      -0.355555_wp, 0.087948_wp, 0.082223_wp, 0.087794_wp, 0.336717_wp, &
                      -0.023286_wp, -0.021755_wp, -0.380101_wp, -0.148638_wp, 0.167518_wp, &
                      0.167134_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.065157_wp, 0.028645_wp, 0.032016_wp, 0.028648_wp, -0.030671_wp, &
                      0.03068_wp, 0.030672_wp, 0.010286_wp, -0.160092_wp, 0.047487_wp, &
                      0.047487_wp]
         case ("resp")
            charges = [ &
                      -0.355208_wp, 0.089729_wp, 0.083442_wp, 0.089639_wp, 0.334726_wp, &
                      -0.023609_wp, -0.022134_wp, -0.333955_wp, -0.205318_wp, 0.171518_wp, &
                      0.17117_wp]
         case ("mbis")
            charges = [ &
                      -0.42907_wp, 0.139466_wp, 0.142676_wp, 0.139465_wp, -0.157985_wp, &
                      0.130392_wp, 0.130402_wp, -0.270512_wp, -0.185828_wp, 0.180496_wp, &
                      0.1805_wp]
         case ("chelpg")
            charges = [ &
                      -0.372368_wp, 0.093716_wp, 0.087202_wp, 0.093604_wp, 0.349569_wp, &
                      -0.026137_wp, -0.024627_wp, -0.344827_wp, -0.204364_wp, 0.174296_wp, &
                      0.173934_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.064851_wp, 0.029297_wp, 0.032283_wp, 0.0293_wp, -0.030386_wp, &
                      0.03079_wp, 0.030781_wp, 0.012273_wp, -0.168051_wp, 0.049283_wp, &
                      0.049283_wp]
         case ("resp")
            charges = [ &
                      -0.35795_wp, 0.090912_wp, 0.084386_wp, 0.090833_wp, 0.334763_wp, &
                      -0.023885_wp, -0.022421_wp, -0.322129_wp, -0.218719_wp, 0.172275_wp, &
                      0.171933_wp]
         case ("mbis")
            charges = [ &
                      -0.430245_wp, 0.140618_wp, 0.143191_wp, 0.140617_wp, -0.158383_wp, &
                      0.130537_wp, 0.130546_wp, -0.266524_wp, -0.196008_wp, 0.182824_wp, &
                      0.182828_wp]
         case ("chelpg")
            charges = [ &
                      -0.375142_wp, 0.094907_wp, 0.088155_wp, 0.094807_wp, 0.349616_wp, &
                      -0.026415_wp, -0.024916_wp, -0.332987_wp, -0.217767_wp, 0.17505_wp, &
                      0.174693_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (26)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.029935_wp, 0.054039_wp, -0.032475_wp, 0.050837_wp, 0.236424_wp, &
                      -0.223516_wp, -0.20739_wp, -0.062957_wp, 0.032936_wp, 0.032663_wp, &
                      0.034328_wp, 0.02768_wp, 0.027498_wp]
         case ("resp")
            charges = [ &
                      -0.294544_wp, 0.125053_wp, 0.285652_wp, 0.116529_wp, 0.758877_wp, &
                      -0.445456_wp, -0.436383_wp, -0.393735_wp, -0.010288_wp, -0.011472_wp, &
                      0.099045_wp, 0.102729_wp, 0.103993_wp]
         case ("mbis")
            charges = [ &
                      -0.248008_wp, 0.162932_wp, -0.192345_wp, 0.164375_wp, 0.658357_wp, &
                      -0.417063_wp, -0.395987_wp, -0.410237_wp, 0.131412_wp, 0.136001_wp, &
                      0.144385_wp, 0.13332_wp, 0.132859_wp]
         case ("chelpg")
            charges = [ &
                      -0.314653_wp, 0.129505_wp, 0.303381_wp, 0.121172_wp, 0.766529_wp, &
                      -0.44692_wp, -0.4383_wp, -0.412151_wp, -0.013043_wp, -0.013996_wp, &
                      0.102718_wp, 0.10724_wp, 0.108519_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.038255_wp, 0.06829_wp, -0.030633_wp, 0.061716_wp, 0.245938_wp, &
                      -0.25238_wp, -0.232048_wp, -0.061522_wp, 0.034026_wp, 0.031078_wp, &
                      0.034615_wp, 0.031405_wp, 0.031261_wp]
         case ("resp")
            charges = [ &
                      -0.260375_wp, 0.138267_wp, 0.31127_wp, 0.12581_wp, 0.803666_wp, &
                      -0.501206_wp, -0.489802_wp, -0.415411_wp, -0.016871_wp, -0.020185_wp, &
                      0.100963_wp, 0.111129_wp, 0.112747_wp]
         case ("mbis")
            charges = [ &
                      -0.260891_wp, 0.187223_wp, -0.182404_wp, 0.183062_wp, 0.706494_wp, &
                      -0.464688_wp, -0.438341_wp, -0.41819_wp, 0.131344_wp, 0.131162_wp, &
                      0.145159_wp, 0.140235_wp, 0.139837_wp]
         case ("chelpg")
            charges = [ &
                      -0.280601_wp, 0.142733_wp, 0.329376_wp, 0.130465_wp, 0.811301_wp, &
                      -0.502642_wp, -0.491725_wp, -0.43407_wp, -0.019719_wp, -0.022804_wp, &
                      0.104681_wp, 0.115685_wp, 0.117319_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.038709_wp, 0.069078_wp, -0.030535_wp, 0.062322_wp, 0.246383_wp, &
                      -0.253851_wp, -0.233385_wp, -0.061418_wp, 0.034036_wp, 0.030949_wp, &
                      0.034582_wp, 0.031636_wp, 0.031494_wp]
         case ("resp")
            charges = [ &
                      -0.258085_wp, 0.138968_wp, 0.312264_wp, 0.126288_wp, 0.805574_wp, &
                      -0.503989_wp, -0.492579_wp, -0.416308_wp, -0.017235_wp, -0.020651_wp, &
                      0.100933_wp, 0.111591_wp, 0.11323_wp]
         case ("mbis")
            charges = [ &
                      -0.261531_wp, 0.18855_wp, -0.181817_wp, 0.184089_wp, 0.70887_wp, &
                      -0.4671_wp, -0.440589_wp, -0.418515_wp, 0.131249_wp, 0.130818_wp, &
                      0.145088_wp, 0.140643_wp, 0.140247_wp]
         case ("chelpg")
            charges = [ &
                      -0.278309_wp, 0.143434_wp, 0.330388_wp, 0.130942_wp, 0.813204_wp, &
                      -0.505418_wp, -0.494508_wp, -0.434979_wp, -0.020088_wp, -0.023274_wp, &
                      0.104653_wp, 0.11615_wp, 0.117805_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (27)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070196_wp, 0.023583_wp, 0.021548_wp, 0.021644_wp, -0.034093_wp, &
                      -0.036281_wp, 0.019326_wp, 0.019183_wp, -0.036059_wp, 0.018311_wp, &
                      0.01849_wp, -0.035953_wp, 0.018729_wp, 0.018407_wp, -0.03517_wp, &
                      0.018481_wp, 0.018904_wp, -0.034287_wp, 0.019894_wp, 0.018802_wp, &
                      -0.043056_wp, 0.019259_wp, 0.021094_wp, 0.041505_wp, 0.02423_wp, &
                      0.018318_wp, -0.250868_wp, 0.018049_wp, 0.031349_wp, 0.14687_wp]
         case ("resp")
            charges = [ &
                      -0.359454_wp, 0.085016_wp, 0.077478_wp, 0.078226_wp, 0.20127_wp, &
                      0.034141_wp, -0.029659_wp, -0.031242_wp, -0.100474_wp, -0.005311_wp, &
                      -0.006861_wp, 0.077111_wp, 0.007739_wp, 0.008447_wp, 0.098332_wp, &
                      -0.02196_wp, -0.017347_wp, -0.023432_wp, -0.021243_wp, -0.026297_wp, &
                      -0.058353_wp, -0.005343_wp, -0.014332_wp, 0.257937_wp, 0.025126_wp, &
                      0.008403_wp, -0.637311_wp, -0.025863_wp, 0.048653_wp, 0.376603_wp]
         case ("mbis")
            charges = [ &
                      -0.421656_wp, 0.130168_wp, 0.127589_wp, 0.127753_wp, -0.160291_wp, &
                      -0.211473_wp, 0.100233_wp, 0.099938_wp, -0.213421_wp, 0.105734_wp, &
                      0.106041_wp, -0.21076_wp, 0.106047_wp, 0.105422_wp, -0.211868_wp, &
                      0.105375_wp, 0.106104_wp, -0.18312_wp, 0.108317_wp, 0.106367_wp, &
                      -0.323063_wp, 0.102752_wp, 0.105765_wp, 0.09511_wp, 0.134723_wp, &
                      0.113077_wp, -0.613541_wp, 0.056314_wp, 0.093438_wp, 0.412943_wp]
         case ("chelpg")
            charges = [ &
                      -0.37774_wp, 0.089572_wp, 0.081486_wp, 0.082232_wp, 0.211754_wp, &
                      0.039691_wp, -0.031523_wp, -0.033133_wp, -0.114375_wp, -0.006402_wp, &
                      -0.007991_wp, 0.086169_wp, 0.010072_wp, 0.01079_wp, 0.103119_wp, &
                      -0.023945_wp, -0.019254_wp, -0.027855_wp, -0.022404_wp, -0.027425_wp, &
                      -0.062527_wp, -0.004577_wp, -0.013749_wp, 0.265846_wp, 0.026123_wp, &
                      0.009205_wp, -0.639521_wp, -0.027687_wp, 0.047126_wp, 0.376926_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072767_wp, 0.023712_wp, 0.021921_wp, 0.021945_wp, -0.034993_wp, &
                      -0.037151_wp, 0.019816_wp, 0.019777_wp, -0.036955_wp, 0.018901_wp, &
                      0.018952_wp, -0.036707_wp, 0.019076_wp, 0.018982_wp, -0.035701_wp, &
                      0.019191_wp, 0.019311_wp, -0.03381_wp, 0.020131_wp, 0.019908_wp, &
                      -0.042622_wp, 0.021682_wp, 0.02228_wp, 0.040155_wp, 0.021786_wp, &
                      0.024425_wp, -0.277085_wp, 0.023699_wp, 0.030698_wp, 0.161461_wp]
         case ("resp")
            charges = [ &
                      -0.380006_wp, 0.088252_wp, 0.080899_wp, 0.081379_wp, 0.221881_wp, &
                      0.040145_wp, -0.034502_wp, -0.036015_wp, -0.100927_wp, -0.008312_wp, &
                      -0.009761_wp, 0.09212_wp, 0.004817_wp, 0.006449_wp, 0.103762_wp, &
                      -0.026097_wp, -0.022865_wp, -0.025006_wp, -0.024005_wp, -0.029343_wp, &
                      -0.042679_wp, -0.004926_wp, -0.014502_wp, 0.2906_wp, 0.0136_wp, &
                      0.014912_wp, -0.725248_wp, -0.020253_wp, 0.040036_wp, 0.425594_wp]
         case ("mbis")
            charges = [ &
                      -0.434824_wp, 0.132995_wp, 0.130845_wp, 0.130875_wp, -0.160633_wp, &
                      -0.214934_wp, 0.101608_wp, 0.101508_wp, -0.216938_wp, 0.107634_wp, &
                      0.107701_wp, -0.213882_wp, 0.107485_wp, 0.107288_wp, -0.215331_wp, &
                      0.107358_wp, 0.107537_wp, -0.185671_wp, 0.109423_wp, 0.108914_wp, &
                      -0.325289_wp, 0.107155_wp, 0.10809_wp, 0.102569_wp, 0.130186_wp, &
                      0.124543_wp, -0.681122_wp, 0.066301_wp, 0.092693_wp, 0.455933_wp]
         case ("chelpg")
            charges = [ &
                      -0.398704_wp, 0.092918_wp, 0.08499_wp, 0.085467_wp, 0.232514_wp, &
                      0.04645_wp, -0.036407_wp, -0.037948_wp, -0.115775_wp, -0.009591_wp, &
                      -0.011074_wp, 0.102214_wp, 0.007262_wp, 0.008908_wp, 0.108698_wp, &
                      -0.028302_wp, -0.024999_wp, -0.030362_wp, -0.02521_wp, -0.030529_wp, &
                      -0.045552_wp, -0.004039_wp, -0.013779_wp, 0.298055_wp, 0.014295_wp, &
                      0.01544_wp, -0.727384_wp, -0.022_wp, 0.038567_wp, 0.425877_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073136_wp, 0.023703_wp, 0.021969_wp, 0.021968_wp, -0.035116_wp, &
                      -0.037269_wp, 0.01987_wp, 0.019869_wp, -0.037076_wp, 0.018995_wp, &
                      0.018995_wp, -0.0368_wp, 0.01909_wp, 0.019082_wp, -0.035752_wp, &
                      0.019333_wp, 0.019329_wp, -0.033689_wp, 0.020085_wp, 0.020142_wp, &
                      -0.042506_wp, 0.022182_wp, 0.022372_wp, 0.039948_wp, 0.021222_wp, &
                      0.025554_wp, -0.281119_wp, 0.024735_wp, 0.030372_wp, 0.163663_wp]
         case ("resp")
            charges = [ &
                      -0.382881_wp, 0.088645_wp, 0.081361_wp, 0.08177_wp, 0.224985_wp, &
                      0.041138_wp, -0.035282_wp, -0.036738_wp, -0.100865_wp, -0.008784_wp, &
                      -0.010277_wp, 0.094403_wp, 0.004273_wp, 0.006141_wp, 0.104523_wp, &
                      -0.026657_wp, -0.023792_wp, -0.025047_wp, -0.024569_wp, -0.029657_wp, &
                      -0.040087_wp, -0.004638_wp, -0.01473_wp, 0.295545_wp, 0.011397_wp, &
                      0.016288_wp, -0.738907_wp, -0.018965_wp, 0.038266_wp, 0.433141_wp]
         case ("mbis")
            charges = [ &
                      -0.436651_wp, 0.133342_wp, 0.131297_wp, 0.13128_wp, -0.160661_wp, &
                      -0.215394_wp, 0.101771_wp, 0.10174_wp, -0.217411_wp, 0.107916_wp, &
                      0.10789_wp, -0.214295_wp, 0.107617_wp, 0.107579_wp, -0.215765_wp, &
                      0.107705_wp, 0.107657_wp, -0.186001_wp, 0.109415_wp, 0.1094_wp, &
                      -0.325439_wp, 0.108044_wp, 0.108254_wp, 0.103802_wp, 0.12909_wp, &
                      0.126628_wp, -0.691464_wp, 0.06814_wp, 0.092153_wp, 0.462377_wp]
         case ("chelpg")
            charges = [ &
                      -0.401633_wp, 0.093324_wp, 0.085462_wp, 0.085869_wp, 0.235642_wp, &
                      0.047557_wp, -0.037193_wp, -0.038679_wp, -0.11585_wp, -0.010092_wp, &
                      -0.011619_wp, 0.10465_wp, 0.006733_wp, 0.008617_wp, 0.109475_wp, &
                      -0.028894_wp, -0.02596_wp, -0.030546_wp, -0.02578_wp, -0.030851_wp, &
                      -0.042739_wp, -0.003733_wp, -0.013987_wp, 0.302926_wp, 0.012037_wp, &
                      0.01677_wp, -0.741036_wp, -0.020697_wp, 0.036802_wp, 0.433422_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (28)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070141_wp, -0.034023_wp, -0.036217_wp, -0.036006_wp, -0.035127_wp, &
                      -0.034281_wp, -0.043032_wp, 0.041502_wp, 0.023698_wp, 0.02174_wp, &
                      0.021597_wp, 0.019421_wp, 0.019235_wp, 0.01841_wp, 0.018725_wp, &
                      0.018943_wp, 0.018517_wp, 0.019918_wp, 0.018823_wp, 0.019272_wp, &
                      0.021103_wp, 0.024233_wp, 0.01832_wp, 0.031361_wp, 0.01802_wp, &
                      -0.250853_wp, 0.146854_wp]
         case ("resp")
            charges = [ &
                      -0.36072_wp, 0.222383_wp, -0.003425_wp, -0.10033_wp, 0.174888_wp, &
                      -0.021494_wp, -0.090289_wp, 0.26959_wp, 0.083807_wp, 0.079266_wp, &
                      0.078307_wp, -0.035279_wp, -0.0337_wp, -0.000187_wp, -0.000544_wp, &
                      0.013158_wp, 0.007401_wp, -0.035781_wp, -0.038615_wp, -0.000508_wp, &
                      -0.010933_wp, 0.030996_wp, 0.01506_wp, 0.045599_wp, -0.028778_wp, &
                      -0.635784_wp, 0.375911_wp]
         case ("mbis")
            charges = [ &
                      -0.421631_wp, -0.160453_wp, -0.210754_wp, -0.216282_wp, -0.210533_wp, &
                      -0.182999_wp, -0.323284_wp, 0.095216_wp, 0.130324_wp, 0.127899_wp, &
                      0.12761_wp, 0.100377_wp, 0.100052_wp, 0.105728_wp, 0.106341_wp, &
                      0.106618_wp, 0.105885_wp, 0.108171_wp, 0.10622_wp, 0.102782_wp, &
                      0.105778_wp, 0.134729_wp, 0.113024_wp, 0.093442_wp, 0.056305_wp, &
                      -0.613527_wp, 0.412976_wp]
         case ("chelpg")
            charges = [ &
                      -0.381153_wp, 0.236708_wp, -0.001433_wp, -0.117768_wp, 0.194977_wp, &
                      -0.027221_wp, -0.100276_wp, 0.28209_wp, 0.088552_wp, 0.083782_wp, &
                      0.082812_wp, -0.037739_wp, -0.0361_wp, -0.000698_wp, -0.001031_wp, &
                      0.016175_wp, 0.010276_wp, -0.039516_wp, -0.042234_wp, 0.000782_wp, &
                      -0.009896_wp, 0.033007_wp, 0.016688_wp, 0.043204_wp, -0.031596_wp, &
                      -0.638861_wp, 0.376467_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072787_wp, -0.034927_wp, -0.037063_wp, -0.036722_wp, -0.035656_wp, &
                      -0.033789_wp, -0.042573_wp, 0.040118_wp, 0.023753_wp, 0.021983_wp, &
                      0.021951_wp, 0.019896_wp, 0.019845_wp, 0.019054_wp, 0.019133_wp, &
                      0.019417_wp, 0.019311_wp, 0.020154_wp, 0.019957_wp, 0.021757_wp, &
                      0.022306_wp, 0.021735_wp, 0.024589_wp, 0.030675_wp, 0.023763_wp, &
                      -0.277585_wp, 0.161716_wp]
         case ("resp")
            charges = [ &
                      -0.380352_wp, 0.243041_wp, 0.000717_wp, -0.095916_wp, 0.185394_wp, &
                      -0.023248_wp, -0.074132_wp, 0.30319_wp, 0.086652_wp, 0.082181_wp, &
                      0.081397_wp, -0.040258_wp, -0.038635_wp, -0.002343_wp, -0.00364_wp, &
                      0.009539_wp, 0.005202_wp, -0.039542_wp, -0.042544_wp, 0.000124_wp, &
                      -0.01107_wp, 0.018914_wp, 0.021533_wp, 0.036564_wp, -0.023337_wp, &
                      -0.725307_wp, 0.425875_wp]
         case ("mbis")
            charges = [ &
                      -0.435031_wp, -0.160752_wp, -0.214247_wp, -0.219644_wp, -0.214024_wp, &
                      -0.185603_wp, -0.325532_wp, 0.102827_wp, 0.133086_wp, 0.130985_wp, &
                      0.130902_wp, 0.101724_wp, 0.101652_wp, 0.107726_wp, 0.107896_wp, &
                      0.108176_wp, 0.108024_wp, 0.109272_wp, 0.108815_wp, 0.107294_wp, &
                      0.108134_wp, 0.130078_wp, 0.124788_wp, 0.092647_wp, 0.06646_wp, &
                      -0.682412_wp, 0.456771_wp]
         case ("chelpg")
            charges = [ &
                      -0.400978_wp, 0.257582_wp, 0.002863_wp, -0.11335_wp, 0.205896_wp, &
                      -0.029602_wp, -0.083363_wp, 0.315505_wp, 0.091438_wp, 0.086738_wp, &
                      0.085942_wp, -0.042783_wp, -0.041102_wp, -0.00293_wp, -0.004207_wp, &
                      0.01252_wp, 0.008048_wp, -0.043361_wp, -0.046259_wp, 0.001503_wp, &
                      -0.009928_wp, 0.020734_wp, 0.022994_wp, 0.034179_wp, -0.026134_wp, &
                      -0.72835_wp, 0.426404_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073117_wp, -0.035033_wp, -0.037161_wp, -0.036793_wp, -0.035697_wp, &
                      -0.033682_wp, -0.042469_wp, 0.039935_wp, 0.023731_wp, 0.02199_wp, &
                      0.02199_wp, 0.019937_wp, 0.019931_wp, 0.019147_wp, 0.01915_wp, &
                      0.019444_wp, 0.019448_wp, 0.020112_wp, 0.020163_wp, 0.022197_wp, &
                      0.022387_wp, 0.021241_wp, 0.025582_wp, 0.030388_wp, 0.024666_wp, &
                      -0.281112_wp, 0.16364_wp]
         case ("resp")
            charges = [ &
                      -0.382713_wp, 0.245727_wp, 0.001324_wp, -0.09526_wp, 0.186673_wp, &
                      -0.023329_wp, -0.071812_wp, 0.307542_wp, 0.086932_wp, 0.082476_wp, &
                      0.081754_wp, -0.040954_wp, -0.039271_wp, -0.002626_wp, -0.004129_wp, &
                      0.008985_wp, 0.00497_wp, -0.040152_wp, -0.042918_wp, 0.000408_wp, &
                      -0.011261_wp, 0.01694_wp, 0.022709_wp, 0.034988_wp, -0.022247_wp, &
                      -0.737232_wp, 0.432475_wp]
         case ("mbis")
            charges = [ &
                      -0.436631_wp, -0.160769_wp, -0.214647_wp, -0.220024_wp, -0.214399_wp, &
                      -0.185891_wp, -0.325661_wp, 0.103905_wp, 0.133369_wp, 0.131319_wp, &
                      0.131297_wp, 0.101854_wp, 0.101864_wp, 0.10799_wp, 0.10802_wp, &
                      0.108297_wp, 0.108348_wp, 0.109261_wp, 0.10924_wp, 0.108075_wp, &
                      0.108276_wp, 0.129117_wp, 0.126621_wp, 0.092171_wp, 0.068064_wp, &
                      -0.691452_wp, 0.462399_wp]
         case ("chelpg")
            charges = [ &
                      -0.403353_wp, 0.260294_wp, 0.003495_wp, -0.112684_wp, 0.207215_wp, &
                      -0.029764_wp, -0.080928_wp, 0.31983_wp, 0.09172_wp, 0.087035_wp, &
                      0.086301_wp, -0.043488_wp, -0.041747_wp, -0.003224_wp, -0.004708_wp, &
                      0.011959_wp, 0.007811_wp, -0.04398_wp, -0.046643_wp, 0.001799_wp, &
                      -0.010106_wp, 0.018729_wp, 0.024146_wp, 0.0326_wp, -0.025038_wp, &
                      -0.740276_wp, 0.433005_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (29)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.034472_wp, 0.019311_wp, 0.021141_wp, -0.032978_wp, 0.020941_wp, &
                      0.019842_wp, -0.069712_wp, 0.022012_wp, 0.022459_wp, 0.024513_wp, &
                      -0.043061_wp, 0.018391_wp, 0.024338_wp, 0.041594_wp, 0.031385_wp, &
                      0.01817_wp, -0.250832_wp, 0.146958_wp]
         case ("resp")
            charges = [ &
                      -0.014955_wp, 0.00452_wp, -0.006951_wp, 0.255122_wp, -0.039405_wp, &
                      -0.040845_wp, -0.370627_wp, 0.079222_wp, 0.084744_wp, 0.08479_wp, &
                      -0.143754_wp, 0.021303_wp, 0.04199_wp, 0.308483_wp, 0.038822_wp, &
                      -0.038131_wp, -0.644057_wp, 0.379729_wp]
         case ("mbis")
            charges = [ &
                      -0.184772_wp, 0.103446_wp, 0.106428_wp, -0.162465_wp, 0.103326_wp, &
                      0.101405_wp, -0.423541_wp, 0.128287_wp, 0.129044_wp, 0.13186_wp, &
                      -0.327116_wp, 0.113557_wp, 0.135227_wp, 0.096494_wp, 0.093283_wp, &
                      0.056219_wp, -0.613678_wp, 0.412997_wp]
         case ("chelpg")
            charges = [ &
                      -0.015954_wp, 0.00493_wp, -0.006801_wp, 0.268844_wp, -0.041834_wp, &
                      -0.043226_wp, -0.389467_wp, 0.083346_wp, 0.088937_wp, 0.089201_wp, &
                      -0.155239_wp, 0.023306_wp, 0.044323_wp, 0.321041_wp, 0.036396_wp, &
                      -0.040996_wp, -0.647073_wp, 0.380266_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.033804_wp, 0.022084_wp, 0.022499_wp, -0.033438_wp, 0.021258_wp, &
                      0.021151_wp, -0.072261_wp, 0.022567_wp, 0.022641_wp, 0.024401_wp, &
                      -0.042535_wp, 0.024977_wp, 0.021683_wp, 0.040193_wp, 0.030589_wp, &
                      0.024311_wp, -0.278804_wp, 0.162487_wp]
         case ("resp")
            charges = [ &
                      -0.017849_wp, 0.006391_wp, -0.005951_wp, 0.275143_wp, -0.043803_wp, &
                      -0.045064_wp, -0.389148_wp, 0.08283_wp, 0.08686_wp, 0.086896_wp, &
                      -0.133005_wp, 0.029402_wp, 0.030744_wp, 0.346971_wp, 0.028694_wp, &
                      -0.032772_wp, -0.73843_wp, 0.432091_wp]
         case ("mbis")
            charges = [ &
                      -0.187493_wp, 0.108435_wp, 0.109032_wp, -0.162835_wp, 0.104238_wp, &
                      0.103958_wp, -0.437157_wp, 0.131931_wp, 0.132018_wp, 0.134368_wp, &
                      -0.329528_wp, 0.12589_wp, 0.130262_wp, 0.10454_wp, 0.092275_wp, &
                      0.067055_wp, -0.685733_wp, 0.458747_wp]
         case ("chelpg")
            charges = [ &
                      -0.019316_wp, 0.006883_wp, -0.005708_wp, 0.289229_wp, -0.046294_wp, &
                      -0.047515_wp, -0.408206_wp, 0.087006_wp, 0.091101_wp, 0.09134_wp, &
                      -0.144152_wp, 0.031339_wp, 0.032993_wp, 0.359531_wp, 0.026253_wp, &
                      -0.035646_wp, -0.741446_wp, 0.432608_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.033722_wp, 0.022389_wp, 0.022561_wp, -0.03346_wp, 0.021238_wp, &
                      0.021297_wp, -0.072455_wp, 0.022625_wp, 0.022624_wp, 0.024358_wp, &
                      -0.042464_wp, 0.025631_wp, 0.021356_wp, 0.040084_wp, 0.030393_wp, &
                      0.024915_wp, -0.281118_wp, 0.163746_wp]
         case ("resp")
            charges = [ &
                      -0.017985_wp, 0.006678_wp, -0.005983_wp, 0.276737_wp, -0.044244_wp, &
                      -0.045323_wp, -0.390575_wp, 0.083141_wp, 0.086961_wp, 0.086994_wp, &
                      -0.131997_wp, 0.030288_wp, 0.029568_wp, 0.350098_wp, 0.027605_wp, &
                      -0.032079_wp, -0.74631_wp, 0.436426_wp]
         case ("mbis")
            charges = [ &
                      -0.187679_wp, 0.108971_wp, 0.109137_wp, -0.162835_wp, 0.104217_wp, &
                      0.104227_wp, -0.438169_wp, 0.132238_wp, 0.13219_wp, 0.134498_wp, &
                      -0.329617_wp, 0.127095_wp, 0.129623_wp, 0.105265_wp, 0.091945_wp, &
                      0.068125_wp, -0.691663_wp, 0.462433_wp]
         case ("chelpg")
            charges = [ &
                      -0.019488_wp, 0.007178_wp, -0.005733_wp, 0.290847_wp, -0.04674_wp, &
                      -0.047778_wp, -0.409648_wp, 0.087321_wp, 0.091204_wp, 0.09144_wp, &
                      -0.143116_wp, 0.032223_wp, 0.031808_wp, 0.362653_wp, 0.025162_wp, &
                      -0.034949_wp, -0.74933_wp, 0.436946_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (30)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.086781_wp, -0.024186_wp, -0.035924_wp, -0.032494_wp, -0.069413_wp, &
                      0.029523_wp, 0.027867_wp, 0.028961_wp, 0.025427_wp, 0.025505_wp, &
                      0.020745_wp, 0.020842_wp, 0.024872_wp, 0.022225_wp, 0.022836_wp]
         case ("resp")
            charges = [ &
                      -0.491807_wp, -0.137012_wp, 0.011252_wp, 0.168408_wp, -0.361324_wp, &
                      0.188906_wp, 0.194725_wp, 0.123563_wp, 0.030938_wp, 0.044512_wp, &
                      -0.00586_wp, -0.016336_wp, 0.090957_wp, 0.080653_wp, 0.078426_wp]
         case ("mbis")
            charges = [ &
                      -0.43129_wp, -0.07134_wp, -0.205922_wp, -0.164217_wp, -0.424505_wp, &
                      0.170357_wp, 0.166619_wp, 0.121745_wp, 0.11532_wp, 0.120021_wp, &
                      0.109377_wp, 0.102932_wp, 0.132681_wp, 0.128347_wp, 0.12988_wp]
         case ("chelpg")
            charges = [ &
                      -0.493505_wp, -0.136834_wp, 0.00912_wp, 0.17809_wp, -0.37622_wp, &
                      0.189555_wp, 0.19539_wp, 0.123597_wp, 0.03134_wp, 0.044763_wp, &
                      -0.007595_wp, -0.018166_wp, 0.094581_wp, 0.084044_wp, 0.081839_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.092608_wp, -0.0263_wp, -0.0367_wp, -0.03261_wp, -0.070305_wp, &
                      0.030415_wp, 0.02934_wp, 0.031146_wp, 0.02699_wp, 0.026712_wp, &
                      0.020277_wp, 0.022052_wp, 0.025276_wp, 0.022952_wp, 0.023367_wp]
         case ("resp")
            charges = [ &
                      -0.514664_wp, -0.148817_wp, 0.013542_wp, 0.1819_wp, -0.372933_wp, &
                      0.196484_wp, 0.204041_wp, 0.131883_wp, 0.032934_wp, 0.045949_wp, &
                      -0.010002_wp, -0.017494_wp, 0.093119_wp, 0.083307_wp, 0.080753_wp]
         case ("mbis")
            charges = [ &
                      -0.452069_wp, -0.074895_wp, -0.210284_wp, -0.163074_wp, -0.431201_wp, &
                      0.176732_wp, 0.173952_wp, 0.127673_wp, 0.118967_wp, 0.123283_wp, &
                      0.108548_wp, 0.105187_wp, 0.13451_wp, 0.130741_wp, 0.131936_wp]
         case ("chelpg")
            charges = [ &
                      -0.516295_wp, -0.148769_wp, 0.011554_wp, 0.191697_wp, -0.387988_wp, &
                      0.19712_wp, 0.204687_wp, 0.131948_wp, 0.033302_wp, 0.046163_wp, &
                      -0.011768_wp, -0.019359_wp, 0.096786_wp, 0.086727_wp, 0.084197_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.099476_wp, -0.028456_wp, -0.037471_wp, -0.032639_wp, -0.071145_wp, &
                      0.031102_wp, 0.030803_wp, 0.033732_wp, 0.02902_wp, 0.02799_wp, &
                      0.019426_wp, 0.023607_wp, 0.025651_wp, 0.023952_wp, 0.02391_wp]
         case ("resp")
            charges = [ &
                      -0.540769_wp, -0.162094_wp, 0.016495_wp, 0.197437_wp, -0.384919_wp, &
                      0.204608_wp, 0.214298_wp, 0.141515_wp, 0.035643_wp, 0.047274_wp, &
                      -0.01552_wp, -0.018501_wp, 0.095115_wp, 0.086409_wp, 0.083011_wp]
         case ("mbis")
            charges = [ &
                      -0.476066_wp, -0.078218_wp, -0.215278_wp, -0.161437_wp, -0.438126_wp, &
                      0.183513_wp, 0.182013_wp, 0.134449_wp, 0.123522_wp, 0.12678_wp, &
                      0.106921_wp, 0.108009_wp, 0.13628_wp, 0.133631_wp, 0.134013_wp]
         case ("chelpg")
            charges = [ &
                      -0.542338_wp, -0.162187_wp, 0.014682_wp, 0.207354_wp, -0.400143_wp, &
                      0.20523_wp, 0.214925_wp, 0.141616_wp, 0.035971_wp, 0.047442_wp, &
                      -0.017323_wp, -0.020402_wp, 0.098827_wp, 0.08986_wp, 0.086487_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (31)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.250743_wp, 0.04147_wp, -0.040688_wp, -0.067444_wp, 0.023172_wp, &
                      0.025409_wp, 0.019467_wp, 0.147225_wp, 0.018272_wp, 0.031525_wp, &
                      0.025059_wp, 0.027273_wp]
         case ("resp")
            charges = [ &
                      -0.624357_wp, 0.180504_wp, 0.083099_wp, -0.319472_wp, 0.073294_wp, &
                      0.012583_wp, -0.001627_wp, 0.371772_wp, -0.002704_wp, 0.067221_wp, &
                      0.067866_wp, 0.091819_wp]
         case ("mbis")
            charges = [ &
                      -0.614792_wp, 0.090616_wp, -0.268922_wp, -0.3963_wp, 0.12574_wp, &
                      0.128828_wp, 0.107522_wp, 0.413037_wp, 0.057523_wp, 0.094403_wp, &
                      0.12886_wp, 0.133483_wp]
         case ("chelpg")
            charges = [ &
                      -0.625425_wp, 0.183074_wp, 0.088901_wp, -0.331812_wp, 0.076174_wp, &
                      0.011427_wp, -0.002796_wp, 0.371791_wp, -0.003553_wp, 0.066458_wp, &
                      0.070713_wp, 0.095048_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.279261_wp, 0.04013_wp, -0.040009_wp, -0.068337_wp, 0.025637_wp, &
                      0.022685_wp, 0.02632_wp, 0.163163_wp, 0.024632_wp, 0.030726_wp, &
                      0.025956_wp, 0.028356_wp]
         case ("resp")
            charges = [ &
                      -0.720726_wp, 0.214501_wp, 0.109564_wp, -0.347176_wp, 0.080433_wp, &
                      -0.000382_wp, 0.004524_wp, 0.425292_wp, 0.004556_wp, 0.058219_wp, &
                      0.073888_wp, 0.097307_wp]
         case ("mbis")
            charges = [ &
                      -0.68848_wp, 0.098801_wp, -0.268227_wp, -0.408421_wp, 0.131737_wp, &
                      0.123416_wp, 0.119842_wp, 0.459887_wp, 0.068701_wp, 0.093382_wp, &
                      0.132207_wp, 0.137155_wp]
         case ("chelpg")
            charges = [ &
                      -0.721763_wp, 0.216934_wp, 0.116096_wp, -0.360113_wp, 0.083422_wp, &
                      -0.001727_wp, 0.003186_wp, 0.425284_wp, 0.003711_wp, 0.057447_wp, &
                      0.076858_wp, 0.100666_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.280966_wp, 0.040053_wp, -0.039948_wp, -0.068367_wp, 0.025831_wp, &
                      0.022444_wp, 0.026809_wp, 0.1641_wp, 0.025081_wp, 0.030581_wp, &
                      0.025965_wp, 0.028415_wp]
         case ("resp")
            charges = [ &
                      -0.726545_wp, 0.216513_wp, 0.111223_wp, -0.348695_wp, 0.080926_wp, &
                      -0.001342_wp, 0.005046_wp, 0.428499_wp, 0.005162_wp, 0.05749_wp, &
                      0.074132_wp, 0.09759_wp]
         case ("mbis")
            charges = [ &
                      -0.692864_wp, 0.099336_wp, -0.268128_wp, -0.409052_wp, 0.132161_wp, &
                      0.122931_wp, 0.120715_wp, 0.46262_wp, 0.069494_wp, 0.093138_wp, &
                      0.132305_wp, 0.137342_wp]
         case ("chelpg")
            charges = [ &
                      -0.727583_wp, 0.218934_wp, 0.117793_wp, -0.361667_wp, 0.083922_wp, &
                      -0.002698_wp, 0.003702_wp, 0.428493_wp, 0.004321_wp, 0.056717_wp, &
                      0.077109_wp, 0.100956_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (32)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.224927_wp, 0.033199_wp, 0.282492_wp, -0.119707_wp, -0.115547_wp, &
                      -0.114164_wp, 0.162601_wp, 0.042738_wp, 0.053313_wp]
         case ("resp")
            charges = [ &
                      -0.526474_wp, -0.041942_wp, 0.484313_wp, -0.182282_wp, -0.205785_wp, &
                      -0.165218_wp, 0.395105_wp, 0.091134_wp, 0.151149_wp]
         case ("mbis")
            charges = [ &
                      -0.546764_wp, -0.203795_wp, 0.792193_wp, -0.265189_wp, -0.259137_wp, &
                      -0.250578_wp, 0.433824_wp, 0.133929_wp, 0.165516_wp]
         case ("chelpg")
            charges = [ &
                      -0.525565_wp, -0.048453_wp, 0.492151_wp, -0.184_wp, -0.207594_wp, &
                      -0.166916_wp, 0.39528_wp, 0.092523_wp, 0.152574_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.246335_wp, 0.035829_wp, 0.2832_wp, -0.123413_wp, -0.120077_wp, &
                      -0.125854_wp, 0.178711_wp, 0.058125_wp, 0.059812_wp]
         case ("resp")
            charges = [ &
                      -0.613653_wp, 0.016185_wp, 0.486301_wp, -0.191117_wp, -0.215408_wp, &
                      -0.187092_wp, 0.446514_wp, 0.108603_wp, 0.149669_wp]
         case ("mbis")
            charges = [ &
                      -0.606866_wp, -0.202725_wp, 0.801331_wp, -0.271083_wp, -0.266334_wp, &
                      -0.266647_wp, 0.476415_wp, 0.159719_wp, 0.176187_wp]
         case ("chelpg")
            charges = [ &
                      -0.613477_wp, 0.013023_wp, 0.493004_wp, -0.19274_wp, -0.217059_wp, &
                      -0.188676_wp, 0.446679_wp, 0.10906_wp, 0.150185_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.247302_wp, 0.035963_wp, 0.283214_wp, -0.123549_wp, -0.12025_wp, &
                      -0.126423_wp, 0.179407_wp, 0.05888_wp, 0.060058_wp]
         case ("resp")
            charges = [ &
                      -0.617642_wp, 0.019045_wp, 0.486059_wp, -0.191407_wp, -0.215711_wp, &
                      -0.188118_wp, 0.448796_wp, 0.109501_wp, 0.149478_wp]
         case ("mbis")
            charges = [ &
                      -0.609553_wp, -0.202572_wp, 0.801639_wp, -0.271292_wp, -0.266602_wp, &
                      -0.26741_wp, 0.478253_wp, 0.160968_wp, 0.176568_wp]
         case ("chelpg")
            charges = [ &
                      -0.617503_wp, 0.016052_wp, 0.492705_wp, -0.193025_wp, -0.217354_wp, &
                      -0.189696_wp, 0.448961_wp, 0.109911_wp, 0.149948_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (33)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.075044_wp, 0.028936_wp, -0.073182_wp, -0.072752_wp, -0.039493_wp, &
                      -0.002536_wp, -0.071852_wp, -0.072345_wp, 0.020796_wp, 0.020973_wp, &
                      0.021046_wp, 0.022249_wp, 0.022158_wp, 0.02219_wp, 0.021886_wp, &
                      0.020878_wp, 0.02124_wp, 0.01752_wp, 0.017414_wp, 0.018608_wp, &
                      0.020948_wp, 0.022021_wp, 0.021897_wp, 0.021985_wp, 0.022174_wp, &
                      0.022282_wp]
         case ("resp")
            charges = [ &
                      -0.543384_wp, 0.861192_wp, -0.546818_wp, -0.515068_wp, -0.468846_wp, &
                      0.614602_wp, -0.544593_wp, -0.525363_wp, 0.103112_wp, 0.101464_wp, &
                      0.111513_wp, 0.106153_wp, 0.101447_wp, 0.106977_wp, 0.096773_wp, &
                      0.097452_wp, 0.103731_wp, 0.062926_wp, 0.072895_wp, -0.077708_wp, &
                      0.123684_wp, 0.116709_wp, 0.111156_wp, 0.104887_wp, 0.109002_wp, &
                      0.116105_wp]
         case ("mbis")
            charges = [ &
                      -0.459931_wp, 0.208178_wp, -0.468307_wp, -0.482149_wp, -0.27088_wp, &
                      0.044516_wp, -0.469271_wp, -0.446206_wp, 0.131701_wp, 0.129617_wp, &
                      0.140645_wp, 0.13744_wp, 0.134071_wp, 0.138507_wp, 0.138641_wp, &
                      0.134838_wp, 0.138129_wp, 0.114089_wp, 0.110883_wp, 0.09791_wp, &
                      0.135749_wp, 0.13545_wp, 0.133079_wp, 0.128138_wp, 0.132042_wp, &
                      0.133122_wp]
         case ("chelpg")
            charges = [ &
                      -0.569073_wp, 0.892794_wp, -0.57095_wp, -0.531368_wp, -0.514246_wp, &
                      0.643073_wp, -0.56283_wp, -0.550589_wp, 0.108547_wp, 0.107339_wp, &
                      0.119027_wp, 0.111324_wp, 0.10761_wp, 0.11309_wp, 0.101424_wp, &
                      0.101063_wp, 0.107052_wp, 0.072697_wp, 0.08286_wp, -0.080102_wp, &
                      0.128698_wp, 0.119918_wp, 0.115561_wp, 0.109716_wp, 0.114696_wp, &
                      0.122669_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076168_wp, 0.029579_wp, -0.074188_wp, -0.073945_wp, -0.039717_wp, &
                      -0.002234_wp, -0.073186_wp, -0.073403_wp, 0.020932_wp, 0.02117_wp, &
                      0.021298_wp, 0.02229_wp, 0.022545_wp, 0.02283_wp, 0.021882_wp, &
                      0.021186_wp, 0.021486_wp, 0.017962_wp, 0.017874_wp, 0.019105_wp, &
                      0.020953_wp, 0.022257_wp, 0.022031_wp, 0.023013_wp, 0.022164_wp, &
                      0.022282_wp]
         case ("resp")
            charges = [ &
                      -0.564174_wp, 0.904732_wp, -0.566899_wp, -0.536974_wp, -0.49317_wp, &
                      0.644985_wp, -0.565121_wp, -0.545322_wp, 0.106141_wp, 0.104501_wp, &
                      0.114852_wp, 0.108696_wp, 0.104584_wp, 0.110629_wp, 0.099783_wp, &
                      0.100954_wp, 0.107434_wp, 0.066314_wp, 0.076632_wp, -0.081961_wp, &
                      0.1274_wp, 0.120399_wp, 0.1146_wp, 0.109393_wp, 0.112083_wp, &
                      0.119509_wp]
         case ("mbis")
            charges = [ &
                      -0.468259_wp, 0.217278_wp, -0.476796_wp, -0.491172_wp, -0.275757_wp, &
                      0.050335_wp, -0.47805_wp, -0.454271_wp, 0.133455_wp, 0.131449_wp, &
                      0.142488_wp, 0.138967_wp, 0.136181_wp, 0.141061_wp, 0.140208_wp, &
                      0.137012_wp, 0.140242_wp, 0.115628_wp, 0.112487_wp, 0.09828_wp, &
                      0.13737_wp, 0.137533_wp, 0.134919_wp, 0.131525_wp, 0.133409_wp, &
                      0.134476_wp]
         case ("chelpg")
            charges = [ &
                      -0.589977_wp, 0.936465_wp, -0.59114_wp, -0.553383_wp, -0.53886_wp, &
                      0.673607_wp, -0.583442_wp, -0.570692_wp, 0.111603_wp, 0.110405_wp, &
                      0.122401_wp, 0.113893_wp, 0.110781_wp, 0.116773_wp, 0.104468_wp, &
                      0.104593_wp, 0.110782_wp, 0.076157_wp, 0.086668_wp, -0.084365_wp, &
                      0.132437_wp, 0.123623_wp, 0.119024_wp, 0.114255_wp, 0.117809_wp, &
                      0.126112_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.077444_wp, 0.030398_wp, -0.07529_wp, -0.075322_wp, -0.039918_wp, &
                      -0.001796_wp, -0.074724_wp, -0.074587_wp, 0.021013_wp, 0.021325_wp, &
                      0.021622_wp, 0.02223_wp, 0.023043_wp, 0.023695_wp, 0.021771_wp, &
                      0.021548_wp, 0.021693_wp, 0.018541_wp, 0.018395_wp, 0.019853_wp, &
                      0.020819_wp, 0.022448_wp, 0.022147_wp, 0.024333_wp, 0.022045_wp, &
                      0.022159_wp]
         case ("resp")
            charges = [ &
                      -0.588224_wp, 0.956905_wp, -0.590188_wp, -0.562775_wp, -0.521866_wp, &
                      0.681602_wp, -0.589131_wp, -0.568543_wp, 0.109412_wp, 0.107771_wp, &
                      0.11865_wp, 0.111338_wp, 0.108182_wp, 0.114986_wp, 0.103077_wp, &
                      0.105039_wp, 0.111584_wp, 0.070309_wp, 0.080905_wp, -0.086874_wp, &
                      0.131401_wp, 0.124523_wp, 0.118537_wp, 0.114854_wp, 0.115396_wp, &
                      0.123131_wp]
         case ("mbis")
            charges = [ &
                      -0.477711_wp, 0.227906_wp, -0.486471_wp, -0.501501_wp, -0.281295_wp, &
                      0.057066_wp, -0.487975_wp, -0.463353_wp, 0.135319_wp, 0.133414_wp, &
                      0.144647_wp, 0.140486_wp, 0.138669_wp, 0.144212_wp, 0.14182_wp, &
                      0.139536_wp, 0.142546_wp, 0.117504_wp, 0.114293_wp, 0.099_wp, &
                      0.138973_wp, 0.139781_wp, 0.136959_wp, 0.135656_wp, 0.134745_wp, &
                      0.135773_wp]
         case ("chelpg")
            charges = [ &
                      -0.614144_wp, 0.988778_wp, -0.614547_wp, -0.579302_wp, -0.56787_wp, &
                      0.710386_wp, -0.607538_wp, -0.594067_wp, 0.114901_wp, 0.113706_wp, &
                      0.126236_wp, 0.116563_wp, 0.114416_wp, 0.121163_wp, 0.107797_wp, &
                      0.108707_wp, 0.114961_wp, 0.080231_wp, 0.09102_wp, -0.089288_wp, &
                      0.136463_wp, 0.127762_wp, 0.122981_wp, 0.119751_wp, 0.121156_wp, &
                      0.129776_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (34)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072066_wp, -0.003158_wp, -0.072839_wp, -0.039013_wp, -0.003168_wp, &
                      -0.072063_wp, -0.072843_wp, 0.020769_wp, 0.022166_wp, 0.022063_wp, &
                      0.018637_wp, 0.021306_wp, 0.022245_wp, 0.022554_wp, 0.01782_wp, &
                      0.017817_wp, 0.018632_wp, 0.02078_wp, 0.022079_wp, 0.022191_wp, &
                      0.022243_wp, 0.022549_wp, 0.021296_wp]
         case ("resp")
            charges = [ &
                      -0.546639_wp, 0.572358_wp, -0.50435_wp, -0.395516_wp, 0.566616_wp, &
                      -0.551719_wp, -0.514652_wp, 0.125886_wp, 0.108908_wp, 0.119567_wp, &
                      -0.062978_wp, 0.108167_wp, 0.099261_wp, 0.108367_wp, 0.069043_wp, &
                      0.070093_wp, -0.057998_wp, 0.128456_wp, 0.121627_wp, 0.109908_wp, &
                      0.10315_wp, 0.110301_wp, 0.112142_wp]
         case ("mbis")
            charges = [ &
                      -0.464037_wp, 0.046631_wp, -0.461763_wp, -0.267719_wp, 0.046549_wp, &
                      -0.464154_wp, -0.46169_wp, 0.133848_wp, 0.133365_wp, 0.134822_wp, &
                      0.096243_wp, 0.133594_wp, 0.134369_wp, 0.135849_wp, 0.110889_wp, &
                      0.110906_wp, 0.096283_wp, 0.133895_wp, 0.134906_wp, 0.133485_wp, &
                      0.134372_wp, 0.13579_wp, 0.133565_wp]
         case ("chelpg")
            charges = [ &
                      -0.563665_wp, 0.596399_wp, -0.529424_wp, -0.426629_wp, 0.590486_wp, &
                      -0.568931_wp, -0.54032_wp, 0.130279_wp, 0.112969_wp, 0.122949_wp, &
                      -0.065345_wp, 0.114369_wp, 0.105247_wp, 0.113856_wp, 0.075552_wp, &
                      0.076604_wp, -0.06023_wp, 0.132931_wp, 0.125089_wp, 0.114017_wp, &
                      0.109331_wp, 0.115974_wp, 0.118492_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073356_wp, -0.002883_wp, -0.073861_wp, -0.039352_wp, -0.002889_wp, &
                      -0.073351_wp, -0.073867_wp, 0.020871_wp, 0.022269_wp, 0.02227_wp, &
                      0.019252_wp, 0.021498_wp, 0.022801_wp, 0.022616_wp, 0.018186_wp, &
                      0.018185_wp, 0.019236_wp, 0.02088_wp, 0.022294_wp, 0.022292_wp, &
                      0.022806_wp, 0.022621_wp, 0.021481_wp]
         case ("resp")
            charges = [ &
                      -0.56689_wp, 0.599419_wp, -0.521626_wp, -0.414629_wp, 0.593888_wp, &
                      -0.571964_wp, -0.532187_wp, 0.129808_wp, 0.112147_wp, 0.123356_wp, &
                      -0.066242_wp, 0.111364_wp, 0.102568_wp, 0.11104_wp, 0.07218_wp, &
                      0.073237_wp, -0.061328_wp, 0.132356_wp, 0.12543_wp, 0.113127_wp, &
                      0.106519_wp, 0.113045_wp, 0.115383_wp]
         case ("mbis")
            charges = [ &
                      -0.472478_wp, 0.051918_wp, -0.469132_wp, -0.272233_wp, 0.051879_wp, &
                      -0.472599_wp, -0.469097_wp, 0.135591_wp, 0.135083_wp, 0.136776_wp, &
                      0.096908_wp, 0.135226_wp, 0.13662_wp, 0.137261_wp, 0.112283_wp, &
                      0.112304_wp, 0.09692_wp, 0.135635_wp, 0.136874_wp, 0.135199_wp, &
                      0.136643_wp, 0.137224_wp, 0.135191_wp]
         case ("chelpg")
            charges = [ &
                      -0.583995_wp, 0.623586_wp, -0.546825_wp, -0.44595_wp, 0.617886_wp, &
                      -0.589253_wp, -0.55798_wp, 0.134221_wp, 0.116227_wp, 0.126753_wp, &
                      -0.068621_wp, 0.117598_wp, 0.108585_wp, 0.116555_wp, 0.078736_wp, &
                      0.079795_wp, -0.063573_wp, 0.13685_wp, 0.128907_wp, 0.117255_wp, &
                      0.112732_wp, 0.118743_wp, 0.121766_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074913_wp, -0.002465_wp, -0.075041_wp, -0.039724_wp, -0.002465_wp, &
                      -0.074906_wp, -0.075049_wp, 0.020895_wp, 0.022314_wp, 0.022449_wp, &
                      0.020137_wp, 0.021675_wp, 0.023592_wp, 0.022582_wp, 0.018616_wp, &
                      0.018619_wp, 0.020107_wp, 0.020901_wp, 0.022486_wp, 0.022332_wp, &
                      0.023606_wp, 0.0226_wp, 0.02165_wp]
         case ("resp")
            charges = [ &
                      -0.591643_wp, 0.633641_wp, -0.54267_wp, -0.438262_wp, 0.628384_wp, &
                      -0.596705_wp, -0.553568_wp, 0.134343_wp, 0.115926_wp, 0.127815_wp, &
                      -0.070253_wp, 0.115073_wp, 0.106754_wp, 0.114009_wp, 0.075911_wp, &
                      0.076977_wp, -0.06543_wp, 0.136861_wp, 0.129911_wp, 0.116877_wp, &
                      0.110789_wp, 0.116108_wp, 0.11915_wp]
         case ("mbis")
            charges = [ &
                      -0.482473_wp, 0.058335_wp, -0.477772_wp, -0.277555_wp, 0.058351_wp, &
                      -0.482601_wp, -0.477787_wp, 0.137503_wp, 0.137001_wp, 0.138995_wp, &
                      0.097939_wp, 0.137047_wp, 0.139509_wp, 0.138709_wp, 0.113912_wp, &
                      0.11394_wp, 0.097915_wp, 0.137541_wp, 0.139118_wp, 0.137107_wp, &
                      0.139558_wp, 0.138702_wp, 0.137001_wp]
         case ("chelpg")
            charges = [ &
                      -0.608833_wp, 0.657951_wp, -0.568008_wp, -0.469814_wp, 0.652526_wp, &
                      -0.614078_wp, -0.579502_wp, 0.138778_wp, 0.120028_wp, 0.131229_wp, &
                      -0.072647_wp, 0.121343_wp, 0.112808_wp, 0.119552_wp, 0.08252_wp, &
                      0.083588_wp, -0.06769_wp, 0.141377_wp, 0.133405_wp, 0.121025_wp, &
                      0.117038_wp, 0.121835_wp, 0.125569_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (35)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.201908_wp, 0.074544_wp, -0.070281_wp, -0.054754_wp, 0.025993_wp, &
                      -0.064692_wp, -0.054551_wp, 0.034747_wp, 0.028183_wp, 0.036203_wp, &
                      0.036197_wp, 0.032837_wp, 0.035657_wp, 0.03543_wp, 0.033457_wp, &
                      0.034671_wp, 0.03827_wp]
         case ("resp")
            charges = [ &
                      -0.7099_wp, 0.868219_wp, -0.676009_wp, -0.790952_wp, 0.698144_wp, &
                      -0.547369_wp, -0.676368_wp, 0.386262_wp, 0.155404_wp, 0.17585_wp, &
                      0.17547_wp, 0.232155_wp, 0.142714_wp, 0.142331_wp, 0.142631_wp, &
                      0.214213_wp, 0.067206_wp]
         case ("mbis")
            charges = [ &
                      -0.499524_wp, 0.449171_wp, -0.53932_wp, -0.387558_wp, 0.256669_wp, &
                      -0.478791_wp, -0.351798_wp, 0.195912_wp, 0.148911_wp, 0.166901_wp, &
                      0.166754_wp, 0.161829_wp, 0.153789_wp, 0.153502_wp, 0.150897_wp, &
                      0.159224_wp, 0.093436_wp]
         case ("chelpg")
            charges = [ &
                      -0.722127_wp, 0.891641_wp, -0.68962_wp, -0.822826_wp, 0.733672_wp, &
                      -0.572958_wp, -0.703724_wp, 0.403527_wp, 0.158463_wp, 0.178551_wp, &
                      0.178174_wp, 0.238983_wp, 0.148007_wp, 0.147631_wp, 0.147902_wp, &
                      0.219853_wp, 0.064853_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.236913_wp, 0.068269_wp, -0.073495_wp, -0.05353_wp, 0.029195_wp, &
                      -0.063652_wp, -0.055067_wp, 0.025373_wp, 0.036612_wp, 0.035928_wp, &
                      0.036025_wp, 0.043717_wp, 0.041003_wp, 0.040782_wp, 0.038524_wp, &
                      0.04468_wp, 0.042556_wp]
         case ("resp")
            charges = [ &
                      -0.825225_wp, 0.926038_wp, -0.72_wp, -0.828451_wp, 0.733406_wp, &
                      -0.579168_wp, -0.719082_wp, 0.426349_wp, 0.179844_wp, 0.183793_wp, &
                      0.183609_wp, 0.257036_wp, 0.157283_wp, 0.156798_wp, 0.156096_wp, &
                      0.240197_wp, 0.071477_wp]
         case ("mbis")
            charges = [ &
                      -0.571422_wp, 0.467253_wp, -0.564028_wp, -0.401899_wp, 0.270392_wp, &
                      -0.497657_wp, -0.365263_wp, 0.193563_wp, 0.167867_wp, 0.17106_wp, &
                      0.171084_wp, 0.183075_wp, 0.164733_wp, 0.164472_wp, 0.161508_wp, &
                      0.179399_wp, 0.105866_wp]
         case ("chelpg")
            charges = [ &
                      -0.837471_wp, 0.949542_wp, -0.733737_wp, -0.860414_wp, 0.769068_wp, &
                      -0.604985_wp, -0.746487_wp, 0.443619_wp, 0.182936_wp, 0.186523_wp, &
                      0.186342_wp, 0.263884_wp, 0.16263_wp, 0.162152_wp, 0.161422_wp, &
                      0.245847_wp, 0.069129_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.242506_wp, 0.067169_wp, -0.073971_wp, -0.053036_wp, 0.029928_wp, &
                      -0.063352_wp, -0.054905_wp, 0.023738_wp, 0.038061_wp, 0.035662_wp, &
                      0.035778_wp, 0.045498_wp, 0.041858_wp, 0.041639_wp, 0.039295_wp, &
                      0.046231_wp, 0.042917_wp]
         case ("resp")
            charges = [ &
                      -0.843601_wp, 0.934278_wp, -0.7258_wp, -0.832985_wp, 0.738482_wp, &
                      -0.583106_wp, -0.724798_wp, 0.432572_wp, 0.183779_wp, 0.184447_wp, &
                      0.184303_wp, 0.26068_wp, 0.159365_wp, 0.158867_wp, 0.157943_wp, &
                      0.24397_wp, 0.071602_wp]
         case ("mbis")
            charges = [ &
                      -0.582918_wp, 0.46952_wp, -0.567155_wp, -0.403182_wp, 0.27248_wp, &
                      -0.500097_wp, -0.366531_wp, 0.19302_wp, 0.170956_wp, 0.171227_wp, &
                      0.171283_wp, 0.186292_wp, 0.166358_wp, 0.166102_wp, 0.163007_wp, &
                      0.182323_wp, 0.107318_wp]
         case ("chelpg")
            charges = [ &
                      -0.85585_wp, 0.957791_wp, -0.73955_wp, -0.864959_wp, 0.774161_wp, &
                      -0.608952_wp, -0.752211_wp, 0.449846_wp, 0.186875_wp, 0.18718_wp, &
                      0.187039_wp, 0.267532_wp, 0.164719_wp, 0.164228_wp, 0.163276_wp, &
                      0.249622_wp, 0.069253_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (36)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.198902_wp, 0.076004_wp, -0.070456_wp, -0.053585_wp, -0.017026_wp, &
                      -0.053586_wp, 0.076003_wp, -0.070455_wp, 0.035904_wp, 0.035894_wp, &
                      0.028585_wp, 0.034597_wp, 0.042041_wp, 0.034597_wp, 0.035899_wp, &
                      0.035898_wp, 0.028586_wp]
         case ("resp")
            charges = [ &
                      -0.728674_wp, 0.764596_wp, -0.654837_wp, -0.56923_wp, 0.12797_wp, &
                      -0.565492_wp, 0.758621_wp, -0.647941_wp, 0.175055_wp, 0.17494_wp, &
                      0.152411_wp, 0.198738_wp, 0.11798_wp, 0.198203_wp, 0.173561_wp, &
                      0.173616_wp, 0.150482_wp]
         case ("mbis")
            charges = [ &
                      -0.530929_wp, 0.445952_wp, -0.539048_wp, -0.340846_wp, -0.011003_wp, &
                      -0.340844_wp, 0.445947_wp, -0.539049_wp, 0.166573_wp, 0.166252_wp, &
                      0.14968_wp, 0.158485_wp, 0.127836_wp, 0.158485_wp, 0.166456_wp, &
                      0.166368_wp, 0.149682_wp]
         case ("chelpg")
            charges = [ &
                      -0.740663_wp, 0.785496_wp, -0.668755_wp, -0.595355_wp, 0.154014_wp, &
                      -0.591493_wp, 0.779596_wp, -0.662251_wp, 0.178_wp, 0.177866_wp, &
                      0.15552_wp, 0.204265_wp, 0.113143_wp, 0.203674_wp, 0.176591_wp, &
                      0.176646_wp, 0.153706_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.226782_wp, 0.070635_wp, -0.07342_wp, -0.054371_wp, -0.015434_wp, &
                      -0.054371_wp, 0.070633_wp, -0.07342_wp, 0.036003_wp, 0.03615_wp, &
                      0.036246_wp, 0.043889_wp, 0.051949_wp, 0.04389_wp, 0.036054_wp, &
                      0.0361_wp, 0.036247_wp]
         case ("resp")
            charges = [ &
                      -0.839779_wp, 0.825608_wp, -0.695699_wp, -0.615459_wp, 0.152816_wp, &
                      -0.611725_wp, 0.820174_wp, -0.690057_wp, 0.182931_wp, 0.183107_wp, &
                      0.174566_wp, 0.222329_wp, 0.132686_wp, 0.221733_wp, 0.181839_wp, &
                      0.18192_wp, 0.17301_wp]
         case ("mbis")
            charges = [ &
                      -0.593263_wp, 0.462643_wp, -0.562859_wp, -0.356702_wp, -0.012447_wp, &
                      -0.356701_wp, 0.462636_wp, -0.562863_wp, 0.171136_wp, 0.171069_wp, &
                      0.167154_wp, 0.17755_wp, 0.145732_wp, 0.177551_wp, 0.171109_wp, &
                      0.171096_wp, 0.167157_wp]
         case ("chelpg")
            charges = [ &
                      -0.851891_wp, 0.846721_wp, -0.709686_wp, -0.64194_wp, 0.179273_wp, &
                      -0.638081_wp, 0.841365_wp, -0.704446_wp, 0.185886_wp, 0.186043_wp, &
                      0.177693_wp, 0.227929_wp, 0.127761_wp, 0.227276_wp, 0.184881_wp, &
                      0.184961_wp, 0.176255_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.232898_wp, 0.069316_wp, -0.074027_wp, -0.054167_wp, -0.014616_wp, &
                      -0.054168_wp, 0.069314_wp, -0.074028_wp, 0.035746_wp, 0.035931_wp, &
                      0.037998_wp, 0.045904_wp, 0.054111_wp, 0.045905_wp, 0.035811_wp, &
                      0.035867_wp, 0.037998_wp]
         case ("resp")
            charges = [ &
                      -0.864256_wp, 0.838004_wp, -0.702945_wp, -0.624544_wp, 0.159481_wp, &
                      -0.620811_wp, 0.832682_wp, -0.697566_wp, 0.183863_wp, 0.184114_wp, &
                      0.179285_wp, 0.227113_wp, 0.135459_wp, 0.226504_wp, 0.18286_wp, &
                      0.182948_wp, 0.177808_wp]
         case ("mbis")
            charges = [ &
                      -0.606882_wp, 0.465484_wp, -0.566988_wp, -0.359076_wp, -0.011978_wp, &
                      -0.359075_wp, 0.465477_wp, -0.566991_wp, 0.171506_wp, 0.171501_wp, &
                      0.170956_wp, 0.181382_wp, 0.149333_wp, 0.181383_wp, 0.171501_wp, &
                      0.171505_wp, 0.170959_wp]
         case ("chelpg")
            charges = [ &
                      -0.876394_wp, 0.859161_wp, -0.716939_wp, -0.651102_wp, 0.186029_wp, &
                      -0.647244_wp, 0.853917_wp, -0.711964_wp, 0.186818_wp, 0.187049_wp, &
                      0.182414_wp, 0.232729_wp, 0.130514_wp, 0.232062_wp, 0.185903_wp, &
                      0.18599_wp, 0.181056_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (37)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068213_wp, 0.034516_wp, 0.032038_wp, 0.03314_wp, 0.052836_wp, &
                      -0.154203_wp, 0.038437_wp, -0.068226_wp, 0.033115_wp, 0.032051_wp, &
                      0.034512_wp]
         case ("resp")
            charges = [ &
                      -0.43645_wp, 0.12032_wp, 0.128203_wp, 0.140077_wp, 0.264004_wp, &
                      -0.238183_wp, 0.067978_wp, -0.427584_wp, 0.137868_wp, 0.126068_wp, &
                      0.117699_wp]
         case ("mbis")
            charges = [ &
                      -0.465136_wp, 0.165555_wp, 0.146046_wp, 0.163085_wp, 0.081362_wp, &
                      -0.239554_wp, 0.139119_wp, -0.465088_wp, 0.163002_wp, 0.146079_wp, &
                      0.165534_wp]
         case ("chelpg")
            charges = [ &
                      -0.454472_wp, 0.124541_wp, 0.132605_wp, 0.144653_wp, 0.275362_wp, &
                      -0.238863_wp, 0.066818_wp, -0.445081_wp, 0.142311_wp, 0.130378_wp, &
                      0.121748_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.067435_wp, 0.034973_wp, 0.041636_wp, 0.03438_wp, 0.061082_wp, &
                      -0.195973_wp, 0.047812_wp, -0.067455_wp, 0.034352_wp, 0.041653_wp, &
                      0.03498_wp]
         case ("resp")
            charges = [ &
                      -0.469085_wp, 0.125083_wp, 0.1477_wp, 0.14668_wp, 0.335241_wp, &
                      -0.31086_wp, 0.072952_wp, -0.460199_wp, 0.144495_wp, 0.14558_wp, &
                      0.122415_wp]
         case ("mbis")
            charges = [ &
                      -0.482059_wp, 0.167869_wp, 0.164234_wp, 0.166861_wp, 0.109175_wp, &
                      -0.294747_wp, 0.151799_wp, -0.482072_wp, 0.166781_wp, 0.164285_wp, &
                      0.167879_wp]
         case ("chelpg")
            charges = [ &
                      -0.487355_wp, 0.129361_wp, 0.152161_wp, 0.15132_wp, 0.346755_wp, &
                      -0.311547_wp, 0.071775_wp, -0.477943_wp, 0.149002_wp, 0.149949_wp, &
                      0.126522_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.067216_wp, 0.034941_wp, 0.043315_wp, 0.034457_wp, 0.062446_wp, &
                      -0.202709_wp, 0.049294_wp, -0.067237_wp, 0.034429_wp, 0.043333_wp, &
                      0.034949_wp]
         case ("resp")
            charges = [ &
                      -0.473937_wp, 0.125636_wp, 0.150988_wp, 0.147428_wp, 0.347006_wp, &
                      -0.322838_wp, 0.073685_wp, -0.465045_wp, 0.145247_wp, 0.148871_wp, &
                      0.122957_wp]
         case ("mbis")
            charges = [ &
                      -0.484183_wp, 0.16795_wp, 0.167295_wp, 0.167143_wp, 0.113495_wp, &
                      -0.303646_wp, 0.153777_wp, -0.484205_wp, 0.167064_wp, 0.167348_wp, &
                      0.167964_wp]
         case ("chelpg")
            charges = [ &
                      -0.492244_wp, 0.129922_wp, 0.155458_wp, 0.152078_wp, 0.358551_wp, &
                      -0.32353_wp, 0.072503_wp, -0.482825_wp, 0.149764_wp, 0.153249_wp, &
                      0.127073_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (38)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071821_wp, -0.040888_wp, 0.071818_wp, -0.076758_wp, -0.246094_wp, &
                      0.143647_wp, 0.022937_wp, 0.023645_wp, 0.022644_wp, 0.018373_wp, &
                      0.024416_wp, 0.029492_wp, 0.028892_wp, 0.021819_wp, 0.027878_wp]
         case ("resp")
            charges = [ &
                      -0.282689_wp, -0.068567_wp, 0.485533_wp, -0.684586_wp, -0.651202_wp, &
                      0.381772_wp, 0.081108_wp, 0.067142_wp, 0.07597_wp, 0.035532_wp, &
                      0.031915_wp, 0.039945_wp, 0.180948_wp, 0.153378_wp, 0.153801_wp]
         case ("mbis")
            charges = [ &
                      -0.401535_wp, -0.28398_wp, 0.288252_wp, -0.567838_wp, -0.620567_wp, &
                      0.410784_wp, 0.126973_wp, 0.126919_wp, 0.137908_wp, 0.108529_wp, &
                      0.125331_wp, 0.084474_wp, 0.164983_wp, 0.143025_wp, 0.156743_wp]
         case ("chelpg")
            charges = [ &
                      -0.292262_wp, -0.072343_wp, 0.497062_wp, -0.705031_wp, -0.653317_wp, &
                      0.382379_wp, 0.084156_wp, 0.069524_wp, 0.07849_wp, 0.036768_wp, &
                      0.033425_wp, 0.038597_wp, 0.186005_wp, 0.158013_wp, 0.158534_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.075095_wp, -0.039525_wp, 0.071533_wp, -0.077554_wp, -0.272884_wp, &
                      0.158768_wp, 0.025222_wp, 0.023528_wp, 0.01684_wp, 0.025092_wp, &
                      0.029027_wp, 0.029145_wp, 0.025788_wp, 0.027796_wp, 0.03232_wp]
         case ("resp")
            charges = [ &
                      -0.299821_wp, -0.055413_wp, 0.54484_wp, -0.722906_wp, -0.75268_wp, &
                      0.437325_wp, 0.088503_wp, 0.067672_wp, 0.069355_wp, 0.042168_wp, &
                      0.032632_wp, 0.030271_wp, 0.182358_wp, 0.169367_wp, 0.166328_wp]
         case ("mbis")
            charges = [ &
                      -0.412599_wp, -0.290522_wp, 0.310183_wp, -0.584869_wp, -0.692445_wp, &
                      0.455395_wp, 0.133766_wp, 0.129333_wp, 0.130305_wp, 0.120974_wp, &
                      0.133367_wp, 0.082845_wp, 0.161651_wp, 0.155998_wp, 0.16662_wp]
         case ("chelpg")
            charges = [ &
                      -0.309752_wp, -0.058693_wp, 0.556379_wp, -0.74345_wp, -0.754796_wp, &
                      0.437915_wp, 0.091622_wp, 0.070121_wp, 0.071936_wp, 0.043282_wp, &
                      0.034011_wp, 0.028885_wp, 0.187448_wp, 0.174015_wp, 0.171078_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.075359_wp, -0.039356_wp, 0.071488_wp, -0.077572_wp, -0.275028_wp, &
                      0.159937_wp, 0.025455_wp, 0.023439_wp, 0.016282_wp, 0.025741_wp, &
                      0.029462_wp, 0.028984_wp, 0.025431_wp, 0.028376_wp, 0.032721_wp]
         case ("resp")
            charges = [ &
                      -0.301332_wp, -0.053797_wp, 0.549071_wp, -0.725338_wp, -0.76075_wp, &
                      0.441699_wp, 0.089228_wp, 0.067587_wp, 0.068661_wp, 0.042799_wp, &
                      0.032719_wp, 0.029295_wp, 0.182136_wp, 0.170731_wp, 0.16729_wp]
         case ("mbis")
            charges = [ &
                      -0.413352_wp, -0.291006_wp, 0.311943_wp, -0.586028_wp, -0.698146_wp, &
                      0.458839_wp, 0.134387_wp, 0.129365_wp, 0.129504_wp, 0.122147_wp, &
                      0.134099_wp, 0.082478_wp, 0.161142_wp, 0.157182_wp, 0.167448_wp]
         case ("chelpg")
            charges = [ &
                      -0.311303_wp, -0.057017_wp, 0.56061_wp, -0.745889_wp, -0.762869_wp, &
                      0.442291_wp, 0.092356_wp, 0.070043_wp, 0.071247_wp, 0.043899_wp, &
                      0.034083_wp, 0.027901_wp, 0.187229_wp, 0.17538_wp, 0.172041_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (39)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070104_wp, 0.057714_wp, -0.147487_wp, -0.035296_wp, -0.071117_wp, &
                      0.032556_wp, 0.031341_wp, 0.033637_wp, 0.036444_wp, 0.028898_wp, &
                      0.02849_wp, 0.023466_wp, 0.027022_wp, 0.024435_wp]
         case ("resp")
            charges = [ &
                      -0.521945_wp, 0.163602_wp, -0.233477_wp, 0.068973_wp, -0.381656_wp, &
                      0.148904_wp, 0.138702_wp, 0.159339_wp, 0.086117_wp, 0.044503_wp, &
                      0.027898_wp, 0.103491_wp, 0.104638_wp, 0.09091_wp]
         case ("mbis")
            charges = [ &
                      -0.475578_wp, 0.091753_wp, -0.263861_wp, -0.212756_wp, -0.411632_wp, &
                      0.161923_wp, 0.14441_wp, 0.165864_wp, 0.133004_wp, 0.132823_wp, &
                      0.11982_wp, 0.155715_wp, 0.131145_wp, 0.12737_wp]
         case ("chelpg")
            charges = [ &
                      -0.541037_wp, 0.170368_wp, -0.234156_wp, 0.075383_wp, -0.397856_wp, &
                      0.153366_wp, 0.143136_wp, 0.1644_wp, 0.085197_wp, 0.043547_wp, &
                      0.026965_wp, 0.10729_wp, 0.108855_wp, 0.094543_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069835_wp, 0.063098_wp, -0.181446_wp, -0.033783_wp, -0.07309_wp, &
                      0.033406_wp, 0.040206_wp, 0.033715_wp, 0.04512_wp, 0.030367_wp, &
                      0.037583_wp, 0.017803_wp, 0.029597_wp, 0.027257_wp]
         case ("resp")
            charges = [ &
                      -0.553027_wp, 0.211272_wp, -0.295945_wp, 0.083619_wp, -0.400509_wp, &
                      0.154015_wp, 0.158002_wp, 0.163958_wp, 0.09359_wp, 0.041463_wp, &
                      0.03759_wp, 0.096513_wp, 0.112181_wp, 0.097277_wp]
         case ("mbis")
            charges = [ &
                      -0.491471_wp, 0.110856_wp, -0.309753_wp, -0.220963_wp, -0.420859_wp, &
                      0.165064_wp, 0.161629_wp, 0.167636_wp, 0.145939_wp, 0.135879_wp, &
                      0.136344_wp, 0.147857_wp, 0.137652_wp, 0.134189_wp]
         case ("chelpg")
            charges = [ &
                      -0.572392_wp, 0.218157_wp, -0.29665_wp, 0.090408_wp, -0.41706_wp, &
                      0.158531_wp, 0.162492_wp, 0.169093_wp, 0.092622_wp, 0.040405_wp, &
                      0.036564_wp, 0.100379_wp, 0.116476_wp, 0.100975_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069691_wp, 0.064088_wp, -0.187534_wp, -0.033435_wp, -0.073378_wp, &
                      0.033427_wp, 0.041957_wp, 0.033596_wp, 0.04664_wp, 0.030493_wp, &
                      0.039408_wp, 0.01656_wp, 0.030064_wp, 0.027804_wp]
         case ("resp")
            charges = [ &
                      -0.557718_wp, 0.219759_wp, -0.307232_wp, 0.086117_wp, -0.403197_wp, &
                      0.154514_wp, 0.161594_wp, 0.16439_wp, 0.094883_wp, 0.04063_wp, &
                      0.039772_wp, 0.094616_wp, 0.113446_wp, 0.098426_wp]
         case ("mbis")
            charges = [ &
                      -0.493675_wp, 0.114165_wp, -0.31799_wp, -0.22224_wp, -0.422043_wp, &
                      0.165277_wp, 0.164893_wp, 0.167603_wp, 0.148177_wp, 0.136119_wp, &
                      0.139599_wp, 0.145947_wp, 0.138751_wp, 0.135416_wp]
         case ("chelpg")
            charges = [ &
                      -0.577118_wp, 0.226656_wp, -0.307939_wp, 0.092959_wp, -0.419792_wp, &
                      0.159038_wp, 0.166091_wp, 0.169534_wp, 0.093909_wp, 0.039557_wp, &
                      0.038733_wp, 0.09849_wp, 0.117751_wp, 0.102132_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (40)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.035082_wp, 0.025498_wp, 0.024382_wp, -0.040387_wp, 0.032041_wp, &
                      0.035247_wp, 0.179234_wp, -0.273975_wp, -0.072423_wp, 0.035824_wp, &
                      0.041647_wp, 0.0422_wp, -0.036397_wp, 0.01726_wp, 0.01781_wp, &
                      -0.03292_wp, 0.021503_wp, 0.021104_wp, -0.070005_wp, 0.021561_wp, &
                      0.021407_wp, 0.024483_wp]
         case ("resp")
            charges = [ &
                      0.007558_wp, 0.019717_wp, 0.013682_wp, -0.201823_wp, 0.047296_wp, &
                      0.05648_wp, 0.71315_wp, -0.546676_wp, -0.60501_wp, 0.141082_wp, &
                      0.174862_wp, 0.158453_wp, -0.037852_wp, 0.016391_wp, 0.016022_wp, &
                      0.164246_wp, -0.017164_wp, -0.01646_wp, -0.336046_wp, 0.074624_wp, &
                      0.073359_wp, 0.08411_wp]
         case ("mbis")
            charges = [ &
                      -0.181298_wp, 0.123237_wp, 0.117595_wp, -0.410312_wp, 0.147647_wp, &
                      0.148761_wp, 0.653427_wp, -0.511754_wp, -0.634067_wp, 0.172438_wp, &
                      0.192282_wp, 0.179767_wp, -0.217724_wp, 0.104651_wp, 0.105484_wp, &
                      -0.161007_wp, 0.103595_wp, 0.102707_wp, -0.421787_wp, 0.127451_wp, &
                      0.127273_wp, 0.131647_wp]
         case ("chelpg")
            charges = [ &
                      0.020053_wp, 0.017961_wp, 0.011702_wp, -0.225057_wp, 0.052775_wp, &
                      0.061317_wp, 0.729757_wp, -0.549525_wp, -0.62866_wp, 0.14708_wp, &
                      0.180341_wp, 0.163746_wp, -0.04737_wp, 0.018343_wp, 0.018058_wp, &
                      0.175619_wp, -0.019418_wp, -0.018618_wp, -0.34942_wp, 0.077689_wp, &
                      0.076426_wp, 0.087202_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037585_wp, 0.01965_wp, 0.021195_wp, -0.037308_wp, 0.043318_wp, &
                      0.046217_wp, 0.188173_wp, -0.323701_wp, -0.070494_wp, 0.04708_wp, &
                      0.039609_wp, 0.052259_wp, -0.035196_wp, 0.020928_wp, 0.021156_wp, &
                      -0.034013_wp, 0.0211_wp, 0.021077_wp, -0.072379_wp, 0.022324_wp, &
                      0.022291_wp, 0.024311_wp]
         case ("resp")
            charges = [ &
                      0.005049_wp, 0.006683_wp, 0.004848_wp, -0.191318_wp, 0.060267_wp, &
                      0.067776_wp, 0.761815_wp, -0.637212_wp, -0.632136_wp, 0.164671_wp, &
                      0.175465_wp, 0.178423_wp, -0.013074_wp, 0.015061_wp, 0.014402_wp, &
                      0.173535_wp, -0.021605_wp, -0.020118_wp, -0.354413_wp, 0.077888_wp, &
                      0.076917_wp, 0.087076_wp]
         case ("mbis")
            charges = [ &
                      -0.176233_wp, 0.113136_wp, 0.112265_wp, -0.437597_wp, 0.169767_wp, &
                      0.170249_wp, 0.707035_wp, -0.591517_wp, -0.662366_wp, 0.194655_wp, &
                      0.191672_wp, 0.199691_wp, -0.218827_wp, 0.110902_wp, 0.111159_wp, &
                      -0.16262_wp, 0.103607_wp, 0.103383_wp, -0.434832_wp, 0.131279_wp, &
                      0.131321_wp, 0.133883_wp]
         case ("chelpg")
            charges = [ &
                      0.016554_wp, 0.005024_wp, 0.002975_wp, -0.214478_wp, 0.065747_wp, &
                      0.072608_wp, 0.778622_wp, -0.640028_wp, -0.656239_wp, 0.170776_wp, &
                      0.181043_wp, 0.183822_wp, -0.020386_wp, 0.016471_wp, 0.015903_wp, &
                      0.184323_wp, -0.023821_wp, -0.02224_wp, -0.367863_wp, 0.080948_wp, &
                      0.079982_wp, 0.090258_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037905_wp, 0.018712_wp, 0.020638_wp, -0.036876_wp, 0.04486_wp, &
                      0.047652_wp, 0.189112_wp, -0.329723_wp, -0.070203_wp, 0.048593_wp, &
                      0.039166_wp, 0.053515_wp, -0.034985_wp, 0.021522_wp, 0.021669_wp, &
                      -0.034162_wp, 0.020923_wp, 0.020979_wp, -0.072593_wp, 0.022447_wp, &
                      0.022447_wp, 0.024223_wp]
         case ("resp")
            charges = [ &
                      0.003629_wp, 0.00488_wp, 0.003693_wp, -0.188633_wp, 0.061983_wp, &
                      0.069111_wp, 0.766209_wp, -0.647888_wp, -0.634409_wp, 0.16766_wp, &
                      0.175003_wp, 0.180757_wp, -0.00932_wp, 0.015072_wp, 0.01431_wp, &
                      0.174_wp, -0.022233_wp, -0.020581_wp, -0.356001_wp, 0.078224_wp, &
                      0.077329_wp, 0.087207_wp]
         case ("mbis")
            charges = [ &
                      -0.175276_wp, 0.111449_wp, 0.111277_wp, -0.440798_wp, 0.172705_wp, &
                      0.172988_wp, 0.713041_wp, -0.601114_wp, -0.665331_wp, 0.197525_wp, &
                      0.191172_wp, 0.202088_wp, -0.219_wp, 0.111908_wp, 0.112017_wp, &
                      -0.162763_wp, 0.103375_wp, 0.103294_wp, -0.436067_wp, 0.131729_wp, &
                      0.131831_wp, 0.133962_wp]
         case ("chelpg")
            charges = [ &
                      0.014799_wp, 0.00327_wp, 0.001872_wp, -0.211646_wp, 0.067439_wp, &
                      0.073919_wp, 0.783001_wp, -0.65069_wp, -0.658557_wp, 0.173774_wp, &
                      0.180592_wp, 0.186167_wp, -0.016191_wp, 0.016382_wp, 0.015713_wp, &
                      0.184671_wp, -0.024436_wp, -0.022691_wp, -0.369473_wp, 0.081285_wp, &
                      0.080394_wp, 0.090407_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (41)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.273724_wp, 0.17927_wp, -0.04045_wp, -0.035213_wp, -0.034058_wp, &
                      -0.068754_wp, -0.072404_wp, 0.032045_wp, 0.035204_wp, 0.025547_wp, &
                      0.024383_wp, 0.018223_wp, 0.018744_wp, 0.024038_wp, 0.023974_wp, &
                      0.023572_wp, 0.04158_wp, 0.036109_wp, 0.041913_wp]
         case ("resp")
            charges = [ &
                      -0.558498_wp, 0.769695_wp, -0.390334_wp, 0.096514_wp, 0.154946_wp, &
                      -0.398423_wp, -0.599248_wp, 0.082914_wp, 0.091977_wp, 0.015828_wp, &
                      0.010141_wp, -0.013375_wp, -0.013858_wp, 0.099828_wp, 0.092186_wp, &
                      0.091301_wp, 0.168462_wp, 0.143437_wp, 0.156509_wp]
         case ("mbis")
            charges = [ &
                      -0.511869_wp, 0.654772_wp, -0.41561_wp, -0.178589_wp, -0.1669_wp, &
                      -0.421259_wp, -0.633525_wp, 0.147866_wp, 0.149233_wp, 0.123286_wp, &
                      0.118017_wp, 0.098919_wp, 0.099884_wp, 0.130427_wp, 0.131017_wp, &
                      0.13033_wp, 0.192181_wp, 0.172653_wp, 0.179168_wp]
         case ("chelpg")
            charges = [ &
                      -0.562007_wp, 0.791198_wp, -0.423095_wp, 0.1119_wp, 0.162934_wp, &
                      -0.418346_wp, -0.626101_wp, 0.089983_wp, 0.098191_wp, 0.013275_wp, &
                      0.007326_wp, -0.014659_wp, -0.014997_wp, 0.105039_wp, 0.096418_wp, &
                      0.095652_wp, 0.174246_wp, 0.150495_wp, 0.162547_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.324408_wp, 0.188383_wp, -0.037276_wp, -0.037693_wp, -0.032784_wp, &
                      -0.071549_wp, -0.070425_wp, 0.043616_wp, 0.046472_wp, 0.019627_wp, &
                      0.021132_wp, 0.02203_wp, 0.022218_wp, 0.024823_wp, 0.023296_wp, &
                      0.023278_wp, 0.039348_wp, 0.047583_wp, 0.052328_wp]
         case ("resp")
            charges = [ &
                      -0.651249_wp, 0.822292_wp, -0.392117_wp, 0.102812_wp, 0.190525_wp, &
                      -0.431151_wp, -0.62771_wp, 0.098806_wp, 0.106011_wp, 0.001838_wp, &
                      0.000293_wp, -0.015986_wp, -0.016814_wp, 0.106241_wp, 0.096073_wp, &
                      0.096013_wp, 0.168724_wp, 0.167991_wp, 0.177409_wp]
         case ("mbis")
            charges = [ &
                      -0.593179_wp, 0.70942_wp, -0.443622_wp, -0.173236_wp, -0.165012_wp, &
                      -0.436589_wp, -0.662414_wp, 0.170538_wp, 0.171286_wp, 0.113083_wp, &
                      0.112592_wp, 0.105036_wp, 0.105414_wp, 0.134782_wp, 0.132765_wp, &
                      0.132752_wp, 0.191305_wp, 0.195301_wp, 0.199776_wp]
         case ("chelpg")
            charges = [ &
                      -0.654794_wp, 0.844123_wp, -0.425611_wp, 0.118964_wp, 0.198757_wp, &
                      -0.451433_wp, -0.654919_wp, 0.106007_wp, 0.112343_wp, -0.000912_wp, &
                      -0.002719_wp, -0.017363_wp, -0.018043_wp, 0.111549_wp, 0.100366_wp, &
                      0.100426_wp, 0.174573_wp, 0.175156_wp, 0.18353_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.329385_wp, 0.189159_wp, -0.036917_wp, -0.037951_wp, -0.0326_wp, &
                      -0.071804_wp, -0.070184_wp, 0.044899_wp, 0.047669_wp, 0.018854_wp, &
                      0.020669_wp, 0.022526_wp, 0.022648_wp, 0.024934_wp, 0.023126_wp, &
                      0.023172_wp, 0.038965_wp, 0.048832_wp, 0.053387_wp]
         case ("resp")
            charges = [ &
                      -0.660096_wp, 0.826088_wp, -0.390822_wp, 0.102406_wp, 0.194436_wp, &
                      -0.434254_wp, -0.629612_wp, 0.10043_wp, 0.107309_wp, 0.000278_wp, &
                      -0.000755_wp, -0.016089_wp, -0.017006_wp, 0.106918_wp, 0.096257_wp, &
                      0.096336_wp, 0.1683_wp, 0.170487_wp, 0.17939_wp]
         case ("mbis")
            charges = [ &
                      -0.601113_wp, 0.714384_wp, -0.446283_wp, -0.172434_wp, -0.164899_wp, &
                      -0.437856_wp, -0.664871_wp, 0.17298_wp, 0.173567_wp, 0.1117_wp, &
                      0.111773_wp, 0.105843_wp, 0.106102_wp, 0.135239_wp, 0.132718_wp, &
                      0.13282_wp, 0.190865_wp, 0.197669_wp, 0.201797_wp]
         case ("chelpg")
            charges = [ &
                      -0.663643_wp, 0.847938_wp, -0.424351_wp, 0.118597_wp, 0.202703_wp, &
                      -0.454571_wp, -0.656845_wp, 0.107638_wp, 0.113646_wp, -0.002486_wp, &
                      -0.00378_wp, -0.017476_wp, -0.018245_wp, 0.112233_wp, 0.100557_wp, &
                      0.100755_wp, 0.174154_wp, 0.17766_wp, 0.185516_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (42)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.040839_wp, 0.035937_wp, -0.251053_wp, 0.025602_wp, 0.037705_wp, &
                      -0.189744_wp, 0.016779_wp, 0.024908_wp, 0.024403_wp, 0.038362_wp, &
                      0.024737_wp, 0.025325_wp, 0.146199_wp]
         case ("resp")
            charges = [ &
                      0.089369_wp, -0.03166_wp, -0.560641_wp, 0.040495_wp, 0.086899_wp, &
                      -0.245467_wp, -0.206427_wp, 0.057711_wp, 0.086963_wp, 0.121446_wp, &
                      0.099722_wp, 0.096637_wp, 0.364954_wp]
         case ("mbis")
            charges = [ &
                      0.027263_wp, -0.02899_wp, -0.611216_wp, 0.080233_wp, 0.105725_wp, &
                      -0.349856_wp, -0.099117_wp, 0.082074_wp, 0.088802_wp, 0.112527_wp, &
                      0.080085_wp, 0.080194_wp, 0.432274_wp]
         case ("chelpg")
            charges = [ &
                      0.093857_wp, -0.03627_wp, -0.561909_wp, 0.039381_wp, 0.086013_wp, &
                      -0.242723_wp, -0.215424_wp, 0.058793_wp, 0.08809_wp, 0.123719_wp, &
                      0.102168_wp, 0.099019_wp, 0.365287_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.038645_wp, 0.037148_wp, -0.276669_wp, 0.028647_wp, 0.042933_wp, &
                      -0.204177_wp, 0.015107_wp, 0.034999_wp, 0.028585_wp, 0.039985_wp, &
                      0.030648_wp, 0.03074_wp, 0.153408_wp]
         case ("resp")
            charges = [ &
                      0.094123_wp, -0.020783_wp, -0.635795_wp, 0.048065_wp, 0.096908_wp, &
                      -0.290671_wp, -0.215797_wp, 0.075604_wp, 0.094969_wp, 0.126408_wp, &
                      0.113319_wp, 0.110942_wp, 0.402709_wp]
         case ("mbis")
            charges = [ &
                      0.024913_wp, -0.031613_wp, -0.668564_wp, 0.086969_wp, 0.116695_wp, &
                      -0.374734_wp, -0.116373_wp, 0.100402_wp, 0.096442_wp, 0.118441_wp, &
                      0.093357_wp, 0.092686_wp, 0.461375_wp]
         case ("chelpg")
            charges = [ &
                      0.09832_wp, -0.024936_wp, -0.637012_wp, 0.047004_wp, 0.096074_wp, &
                      -0.287966_wp, -0.22504_wp, 0.076578_wp, 0.095983_wp, 0.128758_wp, &
                      0.115819_wp, 0.113381_wp, 0.403037_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.038493_wp, 0.037272_wp, -0.27851_wp, 0.028849_wp, 0.043305_wp, &
                      -0.205237_wp, 0.015007_wp, 0.035818_wp, 0.028886_wp, 0.040021_wp, &
                      0.031129_wp, 0.031143_wp, 0.153824_wp]
         case ("resp")
            charges = [ &
                      0.094318_wp, -0.019679_wp, -0.64117_wp, 0.048631_wp, 0.097649_wp, &
                      -0.29416_wp, -0.216278_wp, 0.077025_wp, 0.095476_wp, 0.12658_wp, &
                      0.114388_wp, 0.111965_wp, 0.405254_wp]
         case ("mbis")
            charges = [ &
                      0.02476_wp, -0.031712_wp, -0.672607_wp, 0.087422_wp, 0.117472_wp, &
                      -0.376621_wp, -0.117493_wp, 0.101861_wp, 0.096963_wp, 0.118702_wp, &
                      0.094394_wp, 0.093588_wp, 0.46327_wp]
         case ("chelpg")
            charges = [ &
                      0.09847_wp, -0.023779_wp, -0.642381_wp, 0.047581_wp, 0.096825_wp, &
                      -0.291464_wp, -0.225547_wp, 0.07799_wp, 0.096479_wp, 0.128937_wp, &
                      0.116895_wp, 0.114412_wp, 0.405582_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (43)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07004_wp, -0.009219_wp, -0.073946_wp, 0.040191_wp, -0.243243_wp, &
                      0.024229_wp, 0.023364_wp, 0.026057_wp, 0.017788_wp, 0.023347_wp, &
                      0.02304_wp, 0.022174_wp, 0.030373_wp, 0.018142_wp, 0.147745_wp]
         case ("resp")
            charges = [ &
                      -0.596905_wp, 0.41225_wp, -0.412467_wp, 0.090497_wp, -0.589696_wp, &
                      0.125668_wp, 0.126966_wp, 0.1554_wp, -0.037204_wp, 0.107923_wp, &
                      0.086851_wp, 0.090383_wp, 0.077152_wp, 0.001981_wp, 0.361201_wp]
         case ("mbis")
            charges = [ &
                      -0.441976_wp, -0.053269_wp, -0.432559_wp, 0.063094_wp, -0.60679_wp, &
                      0.135236_wp, 0.130759_wp, 0.139844_wp, 0.092989_wp, 0.132037_wp, &
                      0.141529_wp, 0.129803_wp, 0.096803_wp, 0.058921_wp, 0.413581_wp]
         case ("chelpg")
            charges = [ &
                      -0.615772_wp, 0.424447_wp, -0.430135_wp, 0.092784_wp, -0.590384_wp, &
                      0.130172_wp, 0.131285_wp, 0.16019_wp, -0.038619_wp, 0.112111_wp, &
                      0.090838_wp, 0.094842_wp, 0.076454_wp, 0.00087_wp, 0.360915_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070981_wp, -0.006724_wp, -0.07716_wp, 0.03923_wp, -0.269415_wp, &
                      0.025053_wp, 0.025832_wp, 0.027122_wp, 0.025356_wp, 0.023132_wp, &
                      0.017408_wp, 0.024498_wp, 0.029562_wp, 0.023828_wp, 0.16326_wp]
         case ("resp")
            charges = [ &
                      -0.653606_wp, 0.47567_wp, -0.452531_wp, 0.112878_wp, -0.682479_wp, &
                      0.138145_wp, 0.139172_wp, 0.166284_wp, -0.035104_wp, 0.11348_wp, &
                      0.083981_wp, 0.101546_wp, 0.069857_wp, 0.008644_wp, 0.414062_wp]
         case ("mbis")
            charges = [ &
                      -0.457424_wp, -0.044053_wp, -0.447568_wp, 0.070057_wp, -0.676631_wp, &
                      0.138909_wp, 0.137118_wp, 0.143876_wp, 0.105263_wp, 0.134859_wp, &
                      0.13468_wp, 0.137273_wp, 0.095743_wp, 0.068954_wp, 0.458945_wp]
         case ("chelpg")
            charges = [ &
                      -0.672671_wp, 0.487809_wp, -0.47035_wp, 0.115663_wp, -0.683264_wp, &
                      0.1427_wp, 0.143541_wp, 0.171136_wp, -0.036515_wp, 0.117705_wp, &
                      0.087999_wp, 0.10606_wp, 0.069026_wp, 0.007396_wp, 0.413766_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071011_wp, -0.006503_wp, -0.077402_wp, 0.039162_wp, -0.271372_wp, &
                      0.025062_wp, 0.026069_wp, 0.027209_wp, 0.026032_wp, 0.023033_wp, &
                      0.016908_wp, 0.024712_wp, 0.029379_wp, 0.02433_wp, 0.164393_wp]
         case ("resp")
            charges = [ &
                      -0.657545_wp, 0.480558_wp, -0.455507_wp, 0.114524_wp, -0.689492_wp, &
                      0.138924_wp, 0.140097_wp, 0.167021_wp, -0.034782_wp, 0.113707_wp, &
                      0.083607_wp, 0.102485_wp, 0.069053_wp, 0.009334_wp, 0.418018_wp]
         case ("mbis")
            charges = [ &
                      -0.458393_wp, -0.043433_wp, -0.448566_wp, 0.070675_wp, -0.681821_wp, &
                      0.139049_wp, 0.137649_wp, 0.144149_wp, 0.106377_wp, 0.134908_wp, &
                      0.134012_wp, 0.137898_wp, 0.095424_wp, 0.069838_wp, 0.462236_wp]
         case ("chelpg")
            charges = [ &
                      -0.676623_wp, 0.492691_wp, -0.473333_wp, 0.117338_wp, -0.690283_wp, &
                      0.143482_wp, 0.144469_wp, 0.171876_wp, -0.036192_wp, 0.117934_wp, &
                      0.087626_wp, 0.107001_wp, 0.068215_wp, 0.008078_wp, 0.417721_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (44)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.100723_wp, -0.072756_wp, -0.077951_wp, -0.077953_wp, 0.028911_wp, &
                      0.02767_wp, 0.021942_wp, 0.02193_wp, 0.025598_wp, 0.028929_wp, &
                      0.026735_wp, 0.026715_wp, 0.027679_wp, -0.249877_wp, 0.141713_wp]
         case ("resp")
            charges = [ &
                      0.863553_wp, -0.456498_wp, -0.609787_wp, -0.595883_wp, 0.121622_wp, &
                      0.139026_wp, 0.122595_wp, 0.126002_wp, 0.086146_wp, 0.121001_wp, &
                      0.138279_wp, 0.140522_wp, 0.134043_wp, -0.705235_wp, 0.374613_wp]
         case ("mbis")
            charges = [ &
                      0.482153_wp, -0.50754_wp, -0.564329_wp, -0.564471_wp, 0.153398_wp, &
                      0.162275_wp, 0.141363_wp, 0.141221_wp, 0.143647_wp, 0.153423_wp, &
                      0.157606_wp, 0.157572_wp, 0.16223_wp, -0.630358_wp, 0.41182_wp]
         case ("chelpg")
            charges = [ &
                      0.879122_wp, -0.47369_wp, -0.628872_wp, -0.614764_wp, 0.125712_wp, &
                      0.143571_wp, 0.12704_wp, 0.130507_wp, 0.090546_wp, 0.125129_wp, &
                      0.143146_wp, 0.145451_wp, 0.138563_wp, -0.706337_wp, 0.374875_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.100996_wp, -0.077356_wp, -0.078631_wp, -0.078611_wp, 0.025953_wp, &
                      0.025298_wp, 0.028384_wp, 0.02834_wp, 0.031728_wp, 0.025962_wp, &
                      0.031267_wp, 0.031274_wp, 0.025318_wp, -0.275813_wp, 0.155901_wp]
         case ("resp")
            charges = [ &
                      0.936817_wp, -0.487342_wp, -0.641224_wp, -0.627855_wp, 0.120474_wp, &
                      0.139351_wp, 0.138434_wp, 0.141569_wp, 0.103641_wp, 0.119798_wp, &
                      0.148867_wp, 0.151132_wp, 0.134577_wp, -0.804991_wp, 0.426752_wp]
         case ("mbis")
            charges = [ &
                      0.512839_wp, -0.533856_wp, -0.582254_wp, -0.582312_wp, 0.153357_wp, &
                      0.160164_wp, 0.155099_wp, 0.154918_wp, 0.159545_wp, 0.15337_wp, &
                      0.167571_wp, 0.167601_wp, 0.160117_wp, -0.700919_wp, 0.454768_wp]
         case ("chelpg")
            charges = [ &
                      0.952493_wp, -0.504711_wp, -0.660475_wp, -0.64691_wp, 0.124608_wp, &
                      0.143937_wp, 0.142923_wp, 0.146117_wp, 0.108093_wp, 0.12397_wp, &
                      0.15378_wp, 0.156106_wp, 0.139142_wp, -0.806083_wp, 0.427009_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.101013_wp, -0.077845_wp, -0.078637_wp, -0.078614_wp, 0.025489_wp, &
                      0.024913_wp, 0.029175_wp, 0.029127_wp, 0.032491_wp, 0.025497_wp, &
                      0.031827_wp, 0.031837_wp, 0.024935_wp, -0.278526_wp, 0.157327_wp]
         case ("resp")
            charges = [ &
                      0.943612_wp, -0.490361_wp, -0.643626_wp, -0.63031_wp, 0.120067_wp, &
                      0.138954_wp, 0.140154_wp, 0.143259_wp, 0.105736_wp, 0.119384_wp, &
                      0.149964_wp, 0.152236_wp, 0.134202_wp, -0.815386_wp, 0.432116_wp]
         case ("mbis")
            charges = [ &
                      0.515821_wp, -0.5364_wp, -0.583778_wp, -0.583827_wp, 0.153061_wp, &
                      0.159615_wp, 0.15668_wp, 0.156494_wp, 0.161416_wp, 0.153072_wp, &
                      0.168689_wp, 0.168727_wp, 0.159569_wp, -0.70822_wp, 0.45909_wp]
         case ("chelpg")
            charges = [ &
                      0.959296_wp, -0.507747_wp, -0.662889_wp, -0.649378_wp, 0.124204_wp, &
                      0.143542_wp, 0.144649_wp, 0.147811_wp, 0.110193_wp, 0.12356_wp, &
                      0.154881_wp, 0.157213_wp, 0.138769_wp, -0.81648_wp, 0.432376_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (45)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070047_wp, -0.034883_wp, -0.037595_wp, -0.001496_wp, -0.071804_wp, &
                      -0.072501_wp, 0.02162_wp, 0.023241_wp, 0.021299_wp, 0.018895_wp, &
                      0.019539_wp, 0.017465_wp, 0.018355_wp, 0.017846_wp, 0.021917_wp, &
                      0.022164_wp, 0.020872_wp, 0.020603_wp, 0.022435_wp, 0.022063_wp]
         case ("resp")
            charges = [ &
                      -0.362968_wp, 0.205852_wp, -0.220481_wp, 0.464077_wp, -0.496111_wp, &
                      -0.507216_wp, 0.081214_wp, 0.082067_wp, 0.087855_wp, -0.032714_wp, &
                      -0.0078_wp, 0.050442_wp, 0.045381_wp, -0.049042_wp, 0.098576_wp, &
                      0.118411_wp, 0.108523_wp, 0.121524_wp, 0.117049_wp, 0.09536_wp]
         case ("mbis")
            charges = [ &
                      -0.425098_wp, -0.154663_wp, -0.240893_wp, 0.040639_wp, -0.461708_wp, &
                      -0.449498_wp, 0.128616_wp, 0.130146_wp, 0.127524_wp, 0.096576_wp, &
                      0.103389_wp, 0.108427_wp, 0.109748_wp, 0.090492_wp, 0.131939_wp, &
                      0.135256_wp, 0.13411_wp, 0.130403_wp, 0.133247_wp, 0.131337_wp]
         case ("chelpg")
            charges = [ &
                      -0.377588_wp, 0.225955_wp, -0.243238_wp, 0.485121_wp, -0.513897_wp, &
                      -0.528523_wp, 0.084681_wp, 0.084723_wp, 0.091143_wp, -0.036529_wp, &
                      -0.011741_wp, 0.054895_wp, 0.049668_wp, -0.051344_wp, 0.102605_wp, &
                      0.122195_wp, 0.11306_wp, 0.126657_wp, 0.121651_wp, 0.100506_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07138_wp, -0.035248_wp, -0.037992_wp, -0.001348_wp, -0.073037_wp, &
                      -0.073683_wp, 0.021776_wp, 0.023426_wp, 0.021509_wp, 0.019314_wp, &
                      0.020132_wp, 0.01779_wp, 0.018748_wp, 0.018409_wp, 0.022195_wp, &
                      0.022447_wp, 0.020983_wp, 0.02087_wp, 0.022613_wp, 0.022466_wp]
         case ("resp")
            charges = [ &
                      -0.37426_wp, 0.215191_wp, -0.229355_wp, 0.485425_wp, -0.513842_wp, &
                      -0.524804_wp, 0.083192_wp, 0.084041_wp, 0.090205_wp, -0.034297_wp, &
                      -0.008363_wp, 0.052024_wp, 0.046916_wp, -0.05162_wp, 0.101688_wp, &
                      0.122187_wp, 0.111586_wp, 0.125569_wp, 0.120251_wp, 0.098265_wp]
         case ("mbis")
            charges = [ &
                      -0.432237_wp, -0.15472_wp, -0.244163_wp, 0.045286_wp, -0.469813_wp, &
                      -0.457452_wp, 0.130263_wp, 0.131859_wp, 0.129275_wp, 0.097489_wp, &
                      0.104726_wp, 0.109621_wp, 0.111047_wp, 0.091143_wp, 0.13389_wp, &
                      0.137259_wp, 0.135773_wp, 0.132305_wp, 0.134999_wp, 0.133439_wp]
         case ("chelpg")
            charges = [ &
                      -0.388995_wp, 0.235472_wp, -0.252313_wp, 0.506587_wp, -0.531726_wp, &
                      -0.546237_wp, 0.086684_wp, 0.086717_wp, 0.093517_wp, -0.038145_wp, &
                      -0.012339_wp, 0.056518_wp, 0.051245_wp, -0.053926_wp, 0.105741_wp, &
                      0.125991_wp, 0.116149_wp, 0.130733_wp, 0.124881_wp, 0.103445_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072983_wp, -0.035611_wp, -0.038429_wp, -0.001093_wp, -0.074516_wp, &
                      -0.075077_wp, 0.021874_wp, 0.023578_wp, 0.021704_wp, 0.019905_wp, &
                      0.020902_wp, 0.018149_wp, 0.019228_wp, 0.019219_wp, 0.022509_wp, &
                      0.022734_wp, 0.021021_wp, 0.021153_wp, 0.022736_wp, 0.022988_wp]
         case ("resp")
            charges = [ &
                      -0.387902_wp, 0.2273_wp, -0.240216_wp, 0.51247_wp, -0.535416_wp, &
                      -0.546285_wp, 0.085378_wp, 0.086274_wp, 0.092912_wp, -0.03629_wp, &
                      -0.009163_wp, 0.053788_wp, 0.048689_wp, -0.05479_wp, 0.105401_wp, &
                      0.126642_wp, 0.115051_wp, 0.130396_wp, 0.12392_wp, 0.101843_wp]
         case ("mbis")
            charges = [ &
                      -0.440658_wp, -0.154703_wp, -0.247965_wp, 0.050929_wp, -0.479408_wp, &
                      -0.46685_wp, 0.13207_wp, 0.133778_wp, 0.131266_wp, 0.098705_wp, &
                      0.106397_wp, 0.110964_wp, 0.11259_wp, 0.092135_wp, 0.136183_wp, &
                      0.139559_wp, 0.137581_wp, 0.134499_wp, 0.13691_wp, 0.136006_wp]
         case ("chelpg")
            charges = [ &
                      -0.402758_wp, 0.247783_wp, -0.263388_wp, 0.533776_wp, -0.553398_wp, &
                      -0.567851_wp, 0.088896_wp, 0.08897_wp, 0.096248_wp, -0.040176_wp, &
                      -0.013182_wp, 0.058324_wp, 0.05306_wp, -0.05711_wp, 0.109477_wp, &
                      0.130465_wp, 0.119639_wp, 0.135593_wp, 0.128575_wp, 0.107057_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (46)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.076322_wp, -0.047047_wp, -0.015699_wp, -0.046291_wp, 0.035622_wp, &
                      -0.196336_wp, 0.036287_wp, 0.043569_wp, 0.039527_wp, 0.039844_wp, &
                      -0.068979_wp, 0.037028_wp, 0.037027_wp, 0.029125_wp]
         case ("resp")
            charges = [ &
                      0.848121_wp, -0.588339_wp, 0.194037_wp, -0.493049_wp, 0.371206_wp, &
                      -0.682747_wp, 0.194569_wp, 0.095807_wp, 0.185802_wp, 0.057549_wp, &
                      -0.734035_wp, 0.190465_wp, 0.190109_wp, 0.170505_wp]
         case ("mbis")
            charges = [ &
                      0.436586_wp, -0.31942_wp, -0.021712_wp, -0.284559_wp, 0.178596_wp, &
                      -0.481719_wp, 0.157529_wp, 0.130781_wp, 0.159383_wp, 0.097733_wp, &
                      -0.539738_wp, 0.168084_wp, 0.168117_wp, 0.150337_wp]
         case ("chelpg")
            charges = [ &
                      0.869665_wp, -0.613952_wp, 0.219837_wp, -0.518201_wp, 0.388012_wp, &
                      -0.694095_wp, 0.199898_wp, 0.090873_wp, 0.191354_wp, 0.05505_wp, &
                      -0.74948_wp, 0.193713_wp, 0.193368_wp, 0.173957_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.070959_wp, -0.046199_wp, -0.012511_wp, -0.04641_wp, 0.027307_wp, &
                      -0.229824_wp, 0.047252_wp, 0.054961_wp, 0.049558_wp, 0.044807_wp, &
                      -0.071973_wp, 0.037172_wp, 0.037153_wp, 0.037745_wp]
         case ("resp")
            charges = [ &
                      0.907985_wp, -0.626656_wp, 0.217562_wp, -0.534111_wp, 0.412326_wp, &
                      -0.795768_wp, 0.218651_wp, 0.11239_wp, 0.210734_wp, 0.062655_wp, &
                      -0.779142_wp, 0.199108_wp, 0.198703_wp, 0.195564_wp]
         case ("mbis")
            charges = [ &
                      0.455662_wp, -0.333126_wp, -0.022652_wp, -0.294669_wp, 0.176485_wp, &
                      -0.550732_wp, 0.178679_wp, 0.150625_wp, 0.178612_wp, 0.11088_wp, &
                      -0.565105_wp, 0.172886_wp, 0.172889_wp, 0.169564_wp]
         case ("chelpg")
            charges = [ &
                      0.929649_wp, -0.652405_wp, 0.243501_wp, -0.559398_wp, 0.429211_wp, &
                      -0.807166_wp, 0.224004_wp, 0.10743_wp, 0.216315_wp, 0.060146_wp, &
                      -0.794693_wp, 0.202379_wp, 0.201985_wp, 0.199042_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.070047_wp, -0.045775_wp, -0.011679_wp, -0.046182_wp, 0.025925_wp, &
                      -0.234845_wp, 0.048909_wp, 0.056685_wp, 0.051005_wp, 0.045237_wp, &
                      -0.072391_wp, 0.036984_wp, 0.036962_wp, 0.039116_wp]
         case ("resp")
            charges = [ &
                      0.915937_wp, -0.631447_wp, 0.22189_wp, -0.539615_wp, 0.418243_wp, &
                      -0.812608_wp, 0.221999_wp, 0.114611_wp, 0.214186_wp, 0.062926_wp, &
                      -0.784561_wp, 0.199789_wp, 0.199377_wp, 0.199273_wp]
         case ("mbis")
            charges = [ &
                      0.457901_wp, -0.33436_wp, -0.022186_wp, -0.295522_wp, 0.176007_wp, &
                      -0.56106_wp, 0.181643_wp, 0.153402_wp, 0.181213_wp, 0.112337_wp, &
                      -0.568127_wp, 0.173139_wp, 0.173138_wp, 0.172474_wp]
         case ("chelpg")
            charges = [ &
                      0.937624_wp, -0.657222_wp, 0.247865_wp, -0.56493_wp, 0.435144_wp, &
                      -0.824026_wp, 0.227357_wp, 0.109645_wp, 0.219773_wp, 0.060416_wp, &
                      -0.800121_wp, 0.203063_wp, 0.202661_wp, 0.202753_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (47)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.058221_wp, -0.06649_wp, 0.043738_wp, -0.062293_wp, 0.233563_wp, &
                      -0.220453_wp, -0.212341_wp, 0.037431_wp, 0.037475_wp, 0.035698_wp, &
                      0.038917_wp, 0.039201_wp, 0.037329_wp]
         case ("resp")
            charges = [ &
                      0.282858_wp, -0.468577_wp, 0.019642_wp, -0.436729_wp, 0.650896_wp, &
                      -0.435015_wp, -0.407789_wp, 0.140748_wp, 0.12907_wp, 0.13436_wp, &
                      0.125297_wp, 0.132823_wp, 0.132417_wp]
         case ("mbis")
            charges = [ &
                      0.013523_wp, -0.483815_wp, 0.140083_wp, -0.476133_wp, 0.642321_wp, &
                      -0.417226_wp, -0.398586_wp, 0.159435_wp, 0.1613_wp, 0.1697_wp, &
                      0.163227_wp, 0.161252_wp, 0.164916_wp]
         case ("chelpg")
            charges = [ &
                      0.292225_wp, -0.48755_wp, 0.019149_wp, -0.454106_wp, 0.652789_wp, &
                      -0.436125_wp, -0.408348_wp, 0.145685_wp, 0.133422_wp, 0.139284_wp, &
                      0.129397_wp, 0.137311_wp, 0.136867_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.064948_wp, -0.064958_wp, 0.053383_wp, -0.060139_wp, 0.240905_wp, &
                      -0.247223_wp, -0.236873_wp, 0.048163_wp, 0.04011_wp, 0.033599_wp, &
                      0.040192_wp, 0.049722_wp, 0.038168_wp]
         case ("resp")
            charges = [ &
                      0.344461_wp, -0.495633_wp, 0.026367_wp, -0.472683_wp, 0.685633_wp, &
                      -0.486338_wp, -0.458533_wp, 0.161792_wp, 0.134415_wp, 0.134642_wp, &
                      0.130704_wp, 0.155629_wp, 0.139545_wp]
         case ("mbis")
            charges = [ &
                      0.016761_wp, -0.492616_wp, 0.155926_wp, -0.487346_wp, 0.683702_wp, &
                      -0.460875_wp, -0.439322_wp, 0.178705_wp, 0.165795_wp, 0.166129_wp, &
                      0.165688_wp, 0.180541_wp, 0.166907_wp]
         case ("chelpg")
            charges = [ &
                      0.353874_wp, -0.51479_wp, 0.025891_wp, -0.490277_wp, 0.687523_wp, &
                      -0.487448_wp, -0.459079_wp, 0.166782_wp, 0.13881_wp, 0.139619_wp, &
                      0.134859_wp, 0.160178_wp, 0.144058_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.065277_wp, -0.06486_wp, 0.053884_wp, -0.060012_wp, 0.241202_wp, &
                      -0.248485_wp, -0.2381_wp, 0.048765_wp, 0.040199_wp, 0.033444_wp, &
                      0.040209_wp, 0.050317_wp, 0.038156_wp]
         case ("resp")
            charges = [ &
                      0.347504_wp, -0.496759_wp, 0.026768_wp, -0.474311_wp, 0.68695_wp, &
                      -0.488686_wp, -0.460991_wp, 0.162921_wp, 0.134571_wp, 0.134522_wp, &
                      0.13086_wp, 0.156883_wp, 0.139768_wp]
         case ("mbis")
            charges = [ &
                      0.016865_wp, -0.492896_wp, 0.156752_wp, -0.487757_wp, 0.685571_wp, &
                      -0.462915_wp, -0.441322_wp, 0.179758_wp, 0.165917_wp, 0.165831_wp, &
                      0.165704_wp, 0.181604_wp, 0.166885_wp]
         case ("chelpg")
            charges = [ &
                      0.356921_wp, -0.515924_wp, 0.026294_wp, -0.491913_wp, 0.688836_wp, &
                      -0.489792_wp, -0.461541_wp, 0.167915_wp, 0.138968_wp, 0.139501_wp, &
                      0.135017_wp, 0.161436_wp, 0.144283_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (48)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073236_wp, 0.179299_wp, -0.270559_wp, -0.041412_wp, -0.035753_wp, &
                      -0.036814_wp, -0.035635_wp, -0.034643_wp, -0.070589_wp, 0.041007_wp, &
                      0.038872_wp, 0.03887_wp, 0.033542_wp, 0.033678_wp, 0.024998_wp, &
                      0.024931_wp, 0.017821_wp, 0.017817_wp, 0.020685_wp, 0.020708_wp, &
                      0.019196_wp, 0.019193_wp, 0.022282_wp, 0.023452_wp, 0.022273_wp]
         case ("resp")
            charges = [ &
                      -0.614135_wp, 0.752745_wp, -0.558796_wp, -0.317029_wp, 0.160668_wp, &
                      -0.100044_wp, -0.02983_wp, 0.224895_wp, -0.345066_wp, 0.175673_wp, &
                      0.150619_wp, 0.150473_wp, 0.073109_wp, 0.074973_wp, -0.006895_wp, &
                      -0.004686_wp, 0.018389_wp, 0.017811_wp, 0.010727_wp, 0.008718_wp, &
                      -0.036528_wp, -0.036597_wp, 0.076307_wp, 0.078144_wp, 0.076353_wp]
         case ("mbis")
            charges = [ &
                      -0.630074_wp, 0.657883_wp, -0.511745_wp, -0.409851_wp, -0.177216_wp, &
                      -0.220324_wp, -0.211753_wp, -0.160454_wp, -0.421207_wp, 0.190996_wp, &
                      0.17387_wp, 0.17385_wp, 0.146727_wp, 0.146788_wp, 0.119942_wp, &
                      0.119961_wp, 0.104744_wp, 0.104796_wp, 0.108841_wp, 0.10891_wp, &
                      0.099438_wp, 0.099498_wp, 0.128443_wp, 0.129477_wp, 0.12844_wp]
         case ("chelpg")
            charges = [ &
                      -0.635465_wp, 0.771812_wp, -0.562346_wp, -0.354141_wp, 0.191012_wp, &
                      -0.114019_wp, -0.037578_wp, 0.244853_wp, -0.365859_wp, 0.180236_wp, &
                      0.15592_wp, 0.155802_wp, 0.080683_wp, 0.082496_wp, -0.011857_wp, &
                      -0.00964_wp, 0.021332_wp, 0.020806_wp, 0.011896_wp, 0.009877_wp, &
                      -0.0401_wp, -0.040165_wp, 0.08097_wp, 0.082485_wp, 0.080988_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071542_wp, 0.187865_wp, -0.318043_wp, -0.038451_wp, -0.038289_wp, &
                      -0.035656_wp, -0.03658_wp, -0.035124_wp, -0.073223_wp, 0.038395_wp, &
                      0.04944_wp, 0.04945_wp, 0.044309_wp, 0.044445_wp, 0.02049_wp, &
                      0.020414_wp, 0.02113_wp, 0.021132_wp, 0.020539_wp, 0.020548_wp, &
                      0.02028_wp, 0.020277_wp, 0.022318_wp, 0.023546_wp, 0.022312_wp]
         case ("resp")
            charges = [ &
                      -0.640382_wp, 0.799898_wp, -0.645955_wp, -0.310785_wp, 0.16721_wp, &
                      -0.075289_wp, -0.036557_wp, 0.248632_wp, -0.362377_wp, 0.175162_wp, &
                      0.172206_wp, 0.171998_wp, 0.085492_wp, 0.087431_wp, -0.019607_wp, &
                      -0.017238_wp, 0.015661_wp, 0.015077_wp, 0.008801_wp, 0.00665_wp, &
                      -0.041521_wp, -0.04166_wp, 0.078362_wp, 0.080351_wp, 0.07844_wp]
         case ("mbis")
            charges = [ &
                      -0.658696_wp, 0.710138_wp, -0.588508_wp, -0.436628_wp, -0.172142_wp, &
                      -0.22115_wp, -0.216464_wp, -0.160112_wp, -0.434489_wp, 0.189546_wp, &
                      0.195024_wp, 0.195025_wp, 0.167996_wp, 0.16805_wp, 0.112235_wp, &
                      0.112245_wp, 0.110328_wp, 0.11039_wp, 0.109598_wp, 0.109638_wp, &
                      0.101603_wp, 0.101663_wp, 0.131177_wp, 0.132333_wp, 0.131182_wp]
         case ("chelpg")
            charges = [ &
                      -0.662095_wp, 0.819328_wp, -0.649529_wp, -0.348597_wp, 0.197813_wp, &
                      -0.087865_wp, -0.045533_wp, 0.26927_wp, -0.383362_wp, 0.179787_wp, &
                      0.177611_wp, 0.177431_wp, 0.093185_wp, 0.095072_wp, -0.024685_wp, &
                      -0.022309_wp, 0.01827_wp, 0.017737_wp, 0.01016_wp, 0.007999_wp, &
                      -0.045223_wp, -0.045357_wp, 0.08308_wp, 0.084681_wp, 0.08313_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071215_wp, 0.188995_wp, -0.325285_wp, -0.037924_wp, -0.038701_wp, &
                      -0.035396_wp, -0.036743_wp, -0.035134_wp, -0.073576_wp, 0.037751_wp, &
                      0.051171_wp, 0.051182_wp, 0.046125_wp, 0.046258_wp, 0.019556_wp, &
                      0.019478_wp, 0.021795_wp, 0.021798_wp, 0.020379_wp, 0.020383_wp, &
                      0.020506_wp, 0.020505_wp, 0.022259_wp, 0.023559_wp, 0.022255_wp]
         case ("resp")
            charges = [ &
                      -0.642876_wp, 0.804923_wp, -0.658807_wp, -0.307552_wp, 0.166154_wp, &
                      -0.070229_wp, -0.038184_wp, 0.25219_wp, -0.364556_wp, 0.174358_wp, &
                      0.175477_wp, 0.17526_wp, 0.087359_wp, 0.089303_wp, -0.021607_wp, &
                      -0.019215_wp, 0.015364_wp, 0.014781_wp, 0.008399_wp, 0.006223_wp, &
                      -0.042152_wp, -0.0423_wp, 0.078494_wp, 0.080613_wp, 0.07858_wp]
         case ("mbis")
            charges = [ &
                      -0.662496_wp, 0.717497_wp, -0.600132_wp, -0.440561_wp, -0.170974_wp, &
                      -0.221318_wp, -0.217034_wp, -0.160022_wp, -0.436169_wp, 0.188791_wp, &
                      0.198359_wp, 0.198363_wp, 0.171479_wp, 0.171526_wp, 0.110567_wp, &
                      0.110576_wp, 0.111442_wp, 0.111507_wp, 0.10945_wp, 0.109483_wp, &
                      0.102012_wp, 0.102074_wp, 0.131422_wp, 0.132709_wp, 0.13143_wp]
         case ("chelpg")
            charges = [ &
                      -0.664656_wp, 0.824406_wp, -0.662382_wp, -0.345424_wp, 0.196681_wp, &
                      -0.082428_wp, -0.047443_wp, 0.27294_wp, -0.385542_wp, 0.178994_wp, &
                      0.180899_wp, 0.180711_wp, 0.09506_wp, 0.096952_wp, -0.026683_wp, &
                      -0.024284_wp, 0.017889_wp, 0.017358_wp, 0.009805_wp, 0.007619_wp, &
                      -0.045875_wp, -0.046018_wp, 0.083216_wp, 0.084932_wp, 0.083272_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (49)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069357_wp, -0.032615_wp, -0.040561_wp, 0.179317_wp, -0.072386_wp, &
                      -0.273731_wp, 0.04173_wp, 0.026665_wp, 0.025602_wp, 0.032252_wp, &
                      0.035162_wp, 0.035871_wp, 0.042217_wp, 0.027122_wp, 0.021091_wp, &
                      0.021627_wp]
         case ("resp")
            charges = [ &
                      -0.383894_wp, 0.271571_wp, -0.300299_wp, 0.684053_wp, -0.546611_wp, &
                      -0.544551_wp, 0.163085_wp, -0.016559_wp, -0.020788_wp, 0.070553_wp, &
                      0.08215_wp, 0.128447_wp, 0.146217_wp, 0.089385_wp, 0.089493_wp, &
                      0.087748_wp]
         case ("mbis")
            charges = [ &
                      -0.427106_wp, -0.129699_wp, -0.411472_wp, 0.651074_wp, -0.633255_wp, &
                      -0.511935_wp, 0.19232_wp, 0.11793_wp, 0.112736_wp, 0.147986_wp, &
                      0.149018_wp, 0.172281_wp, 0.179673_wp, 0.135968_wp, 0.126753_wp, &
                      0.127731_wp]
         case ("chelpg")
            charges = [ &
                      -0.402464_wp, 0.294553_wp, -0.330835_wp, 0.702918_wp, -0.570696_wp, &
                      -0.547905_wp, 0.168452_wp, -0.020184_wp, -0.024602_wp, 0.077285_wp, &
                      0.088172_wp, 0.134765_wp, 0.151582_wp, 0.092701_wp, 0.093948_wp, &
                      0.092311_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069791_wp, -0.03516_wp, -0.037297_wp, 0.188512_wp, -0.070387_wp, &
                      -0.3248_wp, 0.039542_wp, 0.02064_wp, 0.022181_wp, 0.043971_wp, &
                      0.046629_wp, 0.047484_wp, 0.052567_wp, 0.026799_wp, 0.024466_wp, &
                      0.024648_wp]
         case ("resp")
            charges = [ &
                      -0.389983_wp, 0.295382_wp, -0.299156_wp, 0.730916_wp, -0.569279_wp, &
                      -0.637102_wp, 0.162467_wp, -0.033501_wp, -0.033597_wp, 0.086464_wp, &
                      0.096526_wp, 0.151741_wp, 0.165583_wp, 0.087267_wp, 0.094285_wp, &
                      0.091988_wp]
         case ("mbis")
            charges = [ &
                      -0.436927_wp, -0.121345_wp, -0.439438_wp, 0.706001_wp, -0.662285_wp, &
                      -0.593881_wp, 0.191524_wp, 0.107285_wp, 0.106713_wp, 0.17083_wp, &
                      0.171355_wp, 0.195192_wp, 0.200169_wp, 0.136541_wp, 0.133961_wp, &
                      0.13431_wp]
         case ("chelpg")
            charges = [ &
                      -0.408853_wp, 0.319098_wp, -0.33033_wp, 0.750076_wp, -0.593701_wp, &
                      -0.640492_wp, 0.167901_wp, -0.037291_wp, -0.037577_wp, 0.093317_wp, &
                      0.102658_wp, 0.158158_wp, 0.171025_wp, 0.090618_wp, 0.09879_wp, &
                      0.096605_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069765_wp, -0.035414_wp, -0.036959_wp, 0.189236_wp, -0.070162_wp, &
                      -0.329434_wp, 0.039192_wp, 0.01992_wp, 0.021741_wp, 0.045167_wp, &
                      0.047753_wp, 0.048656_wp, 0.053537_wp, 0.026675_wp, 0.024869_wp, &
                      0.024991_wp]
         case ("resp")
            charges = [ &
                      -0.389773_wp, 0.296555_wp, -0.297915_wp, 0.733993_wp, -0.570535_wp, &
                      -0.645276_wp, 0.161992_wp, -0.035179_wp, -0.034799_wp, 0.088004_wp, &
                      0.097802_wp, 0.153954_wp, 0.167275_wp, 0.086759_wp, 0.094769_wp, &
                      0.092375_wp]
         case ("mbis")
            charges = [ &
                      -0.437612_wp, -0.120418_wp, -0.441889_wp, 0.710614_wp, -0.664568_wp, &
                      -0.60127_wp, 0.191124_wp, 0.105988_wp, 0.10592_wp, 0.173096_wp, &
                      0.173487_wp, 0.197413_wp, 0.202021_wp, 0.136369_wp, 0.134744_wp, &
                      0.134984_wp]
         case ("chelpg")
            charges = [ &
                      -0.408664_wp, 0.320333_wp, -0.329133_wp, 0.753172_wp, -0.594979_wp, &
                      -0.64867_wp, 0.167429_wp, -0.038985_wp, -0.038794_wp, 0.094866_wp, &
                      0.10394_wp, 0.160379_wp, 0.172721_wp, 0.090112_wp, 0.099278_wp, &
                      0.096995_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (50)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076594_wp, 0.071955_wp, -0.071646_wp, -0.250843_wp, 0.145343_wp, &
                      0.02838_wp, 0.027574_wp, 0.023175_wp, 0.017215_wp, 0.030401_wp, &
                      0.026696_wp, 0.02834_wp]
         case ("resp")
            charges = [ &
                      -0.539878_wp, 0.618366_wp, -0.424635_wp, -0.669086_wp, 0.371809_wp, &
                      0.127797_wp, 0.137274_wp, 0.114662_wp, -0.065032_wp, 0.114869_wp, &
                      0.085373_wp, 0.128481_wp]
         case ("mbis")
            charges = [ &
                      -0.559698_wp, 0.327148_wp, -0.503229_wp, -0.630425_wp, 0.412401_wp, &
                      0.157488_wp, 0.162079_wp, 0.140175_wp, 0.045302_wp, 0.154056_wp, &
                      0.142349_wp, 0.152352_wp]
         case ("chelpg")
            charges = [ &
                      -0.556607_wp, 0.631985_wp, -0.440961_wp, -0.670718_wp, 0.371796_wp, &
                      0.131845_wp, 0.141283_wp, 0.118342_wp, -0.067153_wp, 0.118488_wp, &
                      0.089307_wp, 0.132392_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076702_wp, 0.072289_wp, -0.076338_wp, -0.277951_wp, 0.160698_wp, &
                      0.032932_wp, 0.025484_wp, 0.030259_wp, 0.02385_wp, 0.026597_wp, &
                      0.033101_wp, 0.025778_wp]
         case ("resp")
            charges = [ &
                      -0.558824_wp, 0.671993_wp, -0.470326_wp, -0.767131_wp, 0.425574_wp, &
                      0.135588_wp, 0.135985_wp, 0.128596_wp, -0.058068_wp, 0.11616_wp, &
                      0.107165_wp, 0.133288_wp]
         case ("mbis")
            charges = [ &
                      -0.574189_wp, 0.348853_wp, -0.528953_wp, -0.70348_wp, 0.457983_wp, &
                      0.166632_wp, 0.159722_wp, 0.154217_wp, 0.055205_wp, 0.152431_wp, &
                      0.158612_wp, 0.152965_wp]
         case ("chelpg")
            charges = [ &
                      -0.575671_wp, 0.685646_wp, -0.486907_wp, -0.768739_wp, 0.42555_wp, &
                      0.13967_wp, 0.14003_wp, 0.132307_wp, -0.060169_wp, 0.119843_wp, &
                      0.111173_wp, 0.137267_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076666_wp, 0.072324_wp, -0.076644_wp, -0.279715_wp, 0.161657_wp, &
                      0.033278_wp, 0.02526_wp, 0.030785_wp, 0.024355_wp, 0.026258_wp, &
                      0.033597_wp, 0.025508_wp]
         case ("resp")
            charges = [ &
                      -0.559683_wp, 0.675229_wp, -0.473233_wp, -0.773561_wp, 0.429042_wp, &
                      0.136119_wp, 0.135661_wp, 0.129557_wp, -0.057383_wp, 0.116089_wp, &
                      0.108768_wp, 0.133394_wp]
         case ("mbis")
            charges = [ &
                      -0.574903_wp, 0.350166_wp, -0.530494_wp, -0.708186_wp, 0.460837_wp, &
                      0.167266_wp, 0.159359_wp, 0.155204_wp, 0.055988_wp, 0.152151_wp, &
                      0.159803_wp, 0.152808_wp]
         case ("chelpg")
            charges = [ &
                      -0.576535_wp, 0.688882_wp, -0.489831_wp, -0.775169_wp, 0.42902_wp, &
                      0.140203_wp, 0.139706_wp, 0.133269_wp, -0.059481_wp, 0.119776_wp, &
                      0.112782_wp, 0.137378_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (51)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074631_wp, -0.031848_wp, 0.04179_wp, -0.2476_wp, 0.031122_wp, &
                      0.034834_wp, 0.035903_wp, 0.03799_wp, 0.025163_wp, 0.147278_wp]
         case ("resp")
            charges = [ &
                      -0.476084_wp, -0.177102_wp, 0.267225_wp, -0.630245_wp, 0.197805_wp, &
                      0.197186_wp, 0.152239_wp, 0.06459_wp, 0.007322_wp, 0.397064_wp]
         case ("mbis")
            charges = [ &
                      -0.378117_wp, -0.184144_wp, 0.10777_wp, -0.61985_wp, 0.16217_wp, &
                      0.170441_wp, 0.150792_wp, 0.099433_wp, 0.06769_wp, 0.423816_wp]
         case ("chelpg")
            charges = [ &
                      -0.47742_wp, -0.18048_wp, 0.275559_wp, -0.632622_wp, 0.198404_wp, &
                      0.198055_wp, 0.153047_wp, 0.062726_wp, 0.00514_wp, 0.39759_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07742_wp, -0.035498_wp, 0.040204_wp, -0.273108_wp, 0.037419_wp, &
                      0.038531_wp, 0.037198_wp, 0.039049_wp, 0.032555_wp, 0.161073_wp]
         case ("resp")
            charges = [ &
                      -0.505406_wp, -0.190205_wp, 0.307258_wp, -0.717432_wp, 0.216535_wp, &
                      0.210513_wp, 0.15904_wp, 0.056753_wp, 0.015504_wp, 0.44744_wp]
         case ("mbis")
            charges = [ &
                      -0.39827_wp, -0.192193_wp, 0.112552_wp, -0.68494_wp, 0.177044_wp, &
                      0.180995_wp, 0.156016_wp, 0.102188_wp, 0.081538_wp, 0.465071_wp]
         case ("chelpg")
            charges = [ &
                      -0.506727_wp, -0.193638_wp, 0.315669_wp, -0.719828_wp, 0.217131_wp, &
                      0.211383_wp, 0.159864_wp, 0.05487_wp, 0.013307_wp, 0.447969_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.077464_wp, -0.035767_wp, 0.0401_wp, -0.274757_wp, 0.037885_wp, &
                      0.038769_wp, 0.037196_wp, 0.03902_wp, 0.033088_wp, 0.161932_wp]
         case ("resp")
            charges = [ &
                      -0.506808_wp, -0.191245_wp, 0.310051_wp, -0.723154_wp, 0.217726_wp, &
                      0.211256_wp, 0.159342_wp, 0.056002_wp, 0.016121_wp, 0.450709_wp]
         case ("mbis")
            charges = [ &
                      -0.399147_wp, -0.192852_wp, 0.112956_wp, -0.689108_wp, 0.17801_wp, &
                      0.181578_wp, 0.156218_wp, 0.102183_wp, 0.082523_wp, 0.467641_wp]
         case ("chelpg")
            charges = [ &
                      -0.508128_wp, -0.194681_wp, 0.318468_wp, -0.725551_wp, 0.218322_wp, &
                      0.212126_wp, 0.160167_wp, 0.054117_wp, 0.013923_wp, 0.451238_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (52)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070359_wp, -0.037471_wp, -0.039489_wp, -0.035184_wp, -0.070002_wp, &
                      0.028776_wp, 0.027933_wp, 0.027858_wp, 0.025332_wp, 0.024511_wp, &
                      0.024175_wp, 0.025087_wp, 0.02271_wp, 0.023926_wp, 0.022194_wp]
         case ("resp")
            charges = [ &
                      -0.322218_wp, -0.003496_wp, -0.384033_wp, 0.153282_wp, -0.214036_wp, &
                      0.100775_wp, 0.095718_wp, 0.101745_wp, 0.103055_wp, 0.184511_wp, &
                      0.013779_wp, 0.004477_wp, 0.058241_wp, 0.044181_wp, 0.064018_wp]
         case ("mbis")
            charges = [ &
                      -0.390428_wp, -0.121013_wp, -0.178778_wp, -0.143776_wp, -0.420821_wp, &
                      0.136245_wp, 0.131978_wp, 0.131872_wp, 0.11969_wp, 0.127209_wp, &
                      0.106104_wp, 0.108612_wp, 0.129708_wp, 0.128803_wp, 0.134593_wp]
         case ("chelpg")
            charges = [ &
                      -0.333885_wp, 0.004858_wp, -0.393625_wp, 0.166795_wp, -0.227812_wp, &
                      0.103825_wp, 0.098305_wp, 0.104378_wp, 0.10167_wp, 0.186575_wp, &
                      0.011457_wp, 0.001684_wp, 0.061364_wp, 0.04706_wp, 0.06735_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072063_wp, -0.040828_wp, -0.04292_wp, -0.036221_wp, -0.071082_wp, &
                      0.029732_wp, 0.029156_wp, 0.029124_wp, 0.027168_wp, 0.026296_wp, &
                      0.025229_wp, 0.026152_wp, 0.023675_wp, 0.024975_wp, 0.021601_wp]
         case ("resp")
            charges = [ &
                      -0.331581_wp, -0.011632_wp, -0.405672_wp, 0.161896_wp, -0.216052_wp, &
                      0.104162_wp, 0.099967_wp, 0.106473_wp, 0.110944_wp, 0.194944_wp, &
                      0.014254_wp, 0.004357_wp, 0.059645_wp, 0.045474_wp, 0.062822_wp]
         case ("mbis")
            charges = [ &
                      -0.398937_wp, -0.129094_wp, -0.189821_wp, -0.14507_wp, -0.426882_wp, &
                      0.139723_wp, 0.135835_wp, 0.135804_wp, 0.125717_wp, 0.133599_wp, &
                      0.108684_wp, 0.111303_wp, 0.132647_wp, 0.131774_wp, 0.134713_wp]
         case ("chelpg")
            charges = [ &
                      -0.343035_wp, -0.003556_wp, -0.415091_wp, 0.175473_wp, -0.229907_wp, &
                      0.107169_wp, 0.102516_wp, 0.109069_wp, 0.10963_wp, 0.196968_wp, &
                      0.011902_wp, 0.001541_wp, 0.062778_wp, 0.048375_wp, 0.066169_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073787_wp, -0.044425_wp, -0.046642_wp, -0.037288_wp, -0.072109_wp, &
                      0.030655_wp, 0.030505_wp, 0.030526_wp, 0.029081_wp, 0.028123_wp, &
                      0.026423_wp, 0.02723_wp, 0.024766_wp, 0.026272_wp, 0.020667_wp]
         case ("resp")
            charges = [ &
                      -0.34027_wp, -0.021235_wp, -0.429087_wp, 0.171654_wp, -0.217547_wp, &
                      0.107351_wp, 0.104425_wp, 0.111419_wp, 0.119505_wp, 0.205954_wp, &
                      0.014837_wp, 0.003966_wp, 0.061147_wp, 0.047118_wp, 0.060763_wp]
         case ("mbis")
            charges = [ &
                      -0.407443_wp, -0.137834_wp, -0.201805_wp, -0.146351_wp, -0.432652_wp, &
                      0.143149_wp, 0.139924_wp, 0.139987_wp, 0.132084_wp, 0.140318_wp, &
                      0.111531_wp, 0.114039_wp, 0.135769_wp, 0.135135_wp, 0.134146_wp]
         case ("chelpg")
            charges = [ &
                      -0.351461_wp, -0.013485_wp, -0.438298_wp, 0.185287_wp, -0.231473_wp, &
                      0.110304_wp, 0.106926_wp, 0.113968_wp, 0.118273_wp, 0.207929_wp, &
                      0.012455_wp, 0.001128_wp, 0.064285_wp, 0.050041_wp, 0.064121_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (53)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.004232_wp, 0.027985_wp, -0.192918_wp, 0.030445_wp, -0.040361_wp, &
                      -0.021418_wp, 0.036081_wp, 0.039765_wp, 0.040354_wp, 0.041078_wp, &
                      -0.065595_wp, 0.033203_wp, 0.033199_wp, 0.03395_wp]
         case ("resp")
            charges = [ &
                      0.054568_wp, 0.2536_wp, -0.593491_wp, 0.374079_wp, -0.460817_wp, &
                      0.033111_wp, 0.066767_wp, 0.046916_wp, 0.185609_wp, 0.116976_wp, &
                      -0.535186_wp, 0.155169_wp, 0.154591_wp, 0.148107_wp]
         case ("mbis")
            charges = [ &
                      0.029831_wp, 0.102325_wp, -0.423238_wp, 0.152957_wp, -0.250936_wp, &
                      -0.103316_wp, 0.100102_wp, 0.10026_wp, 0.156082_wp, 0.13791_wp, &
                      -0.442424_wp, 0.147133_wp, 0.147109_wp, 0.146205_wp]
         case ("chelpg")
            charges = [ &
                      0.053153_wp, 0.260257_wp, -0.59943_wp, 0.3846_wp, -0.475149_wp, &
                      0.041721_wp, 0.065601_wp, 0.045012_wp, 0.189233_wp, 0.115949_wp, &
                      -0.551731_wp, 0.15941_wp, 0.158832_wp, 0.152542_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.006267_wp, 0.019404_wp, -0.231674_wp, 0.020922_wp, -0.039656_wp, &
                      -0.016432_wp, 0.041421_wp, 0.044368_wp, 0.051135_wp, 0.053498_wp, &
                      -0.064668_wp, 0.038512_wp, 0.038514_wp, 0.038389_wp]
         case ("resp")
            charges = [ &
                      0.047774_wp, 0.286684_wp, -0.706989_wp, 0.409283_wp, -0.495014_wp, &
                      0.05077_wp, 0.074308_wp, 0.052537_wp, 0.210201_wp, 0.137271_wp, &
                      -0.575275_wp, 0.172699_wp, 0.172114_wp, 0.163636_wp]
         case ("mbis")
            charges = [ &
                      0.038018_wp, 0.0985_wp, -0.497252_wp, 0.150398_wp, -0.259486_wp, &
                      -0.104461_wp, 0.114409_wp, 0.113112_wp, 0.176157_wp, 0.159168_wp, &
                      -0.460526_wp, 0.158058_wp, 0.158047_wp, 0.155858_wp]
         case ("chelpg")
            charges = [ &
                      0.045269_wp, 0.293977_wp, -0.713314_wp, 0.420399_wp, -0.510382_wp, &
                      0.060643_wp, 0.073062_wp, 0.050564_wp, 0.214039_wp, 0.135966_wp, &
                      -0.591401_wp, 0.176872_wp, 0.176289_wp, 0.168016_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.006636_wp, 0.018214_wp, -0.236516_wp, 0.01963_wp, -0.039336_wp, &
                      -0.0155_wp, 0.041838_wp, 0.04469_wp, 0.052467_wp, 0.055116_wp, &
                      -0.064477_wp, 0.039205_wp, 0.039208_wp, 0.038824_wp]
         case ("resp")
            charges = [ &
                      0.046514_wp, 0.290775_wp, -0.721046_wp, 0.413312_wp, -0.498546_wp, &
                      0.053883_wp, 0.074823_wp, 0.052867_wp, 0.213027_wp, 0.139641_wp, &
                      -0.579431_wp, 0.174785_wp, 0.1742_wp, 0.165195_wp]
         case ("mbis")
            charges = [ &
                      0.039016_wp, 0.098033_wp, -0.506541_wp, 0.14995_wp, -0.259978_wp, &
                      -0.103929_wp, 0.11577_wp, 0.114304_wp, 0.178471_wp, 0.16173_wp, &
                      -0.462387_wp, 0.159396_wp, 0.159387_wp, 0.156776_wp]
         case ("chelpg")
            charges = [ &
                      0.04382_wp, 0.298174_wp, -0.727438_wp, 0.424526_wp, -0.514085_wp, &
                      0.063974_wp, 0.073563_wp, 0.050881_wp, 0.2169_wp, 0.138287_wp, &
                      -0.595466_wp, 0.178943_wp, 0.178359_wp, 0.169561_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (54)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.177713_wp, -0.264961_wp, -0.038694_wp, -0.038694_wp, 0.033352_wp, &
                      0.035166_wp, 0.033327_wp, 0.035188_wp, -0.068299_wp, 0.029018_wp, &
                      0.028411_wp, 0.024672_wp, -0.068296_wp, 0.029024_wp, 0.028404_wp, &
                      0.024675_wp]
         case ("resp")
            charges = [ &
                      0.491006_wp, -0.496667_wp, 0.010487_wp, 0.005229_wp, 0.004392_wp, &
                      0.018874_wp, 0.005005_wp, 0.01903_wp, -0.242945_wp, 0.07308_wp, &
                      0.073192_wp, 0.066272_wp, -0.228225_wp, 0.069059_wp, 0.069586_wp, &
                      0.062626_wp]
         case ("mbis")
            charges = [ &
                      0.592461_wp, -0.487219_wp, -0.357542_wp, -0.35755_wp, 0.139704_wp, &
                      0.141515_wp, 0.139668_wp, 0.141553_wp, -0.385643_wp, 0.142989_wp, &
                      0.140242_wp, 0.126118_wp, -0.385636_wp, 0.143016_wp, 0.140211_wp, &
                      0.126121_wp]
         case ("chelpg")
            charges = [ &
                      0.491768_wp, -0.497027_wp, 0.01358_wp, 0.007907_wp, 0.00394_wp, &
                      0.018503_wp, 0.004644_wp, 0.018747_wp, -0.253285_wp, 0.075665_wp, &
                      0.075765_wp, 0.069023_wp, -0.238202_wp, 0.071577_wp, 0.072089_wp, &
                      0.065306_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.187565_wp, -0.313299_wp, -0.035112_wp, -0.035112_wp, 0.045001_wp, &
                      0.046647_wp, 0.044977_wp, 0.046665_wp, -0.071926_wp, 0.023497_wp, &
                      0.024496_wp, 0.030265_wp, -0.071922_wp, 0.023487_wp, 0.024508_wp, &
                      0.030268_wp]
         case ("resp")
            charges = [ &
                      0.537904_wp, -0.589576_wp, 0.020001_wp, 0.016104_wp, 0.017782_wp, &
                      0.032204_wp, 0.017924_wp, 0.031935_wp, -0.253721_wp, 0.063698_wp, &
                      0.066532_wp, 0.079739_wp, -0.239008_wp, 0.059529_wp, 0.062881_wp, &
                      0.076071_wp]
         case ("mbis")
            charges = [ &
                      0.645239_wp, -0.56545_wp, -0.381261_wp, -0.381265_wp, 0.161753_wp, &
                      0.163207_wp, 0.161719_wp, 0.163236_wp, -0.39265_wp, 0.135342_wp, &
                      0.135421_wp, 0.138296_wp, -0.392642_wp, 0.135339_wp, 0.135423_wp, &
                      0.138299_wp]
         case ("chelpg")
            charges = [ &
                      0.538366_wp, -0.589884_wp, 0.023455_wp, 0.019261_wp, 0.017251_wp, &
                      0.031763_wp, 0.017455_wp, 0.031553_wp, -0.264325_wp, 0.066341_wp, &
                      0.069167_wp, 0.082535_wp, -0.249326_wp, 0.062118_wp, 0.065458_wp, &
                      0.078812_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.18828_wp, -0.317311_wp, -0.034771_wp, -0.034771_wp, 0.046063_wp, &
                      0.047668_wp, 0.04604_wp, 0.047685_wp, -0.072203_wp, 0.022899_wp, &
                      0.024052_wp, 0.03081_wp, -0.072198_wp, 0.022887_wp, 0.024066_wp, &
                      0.030812_wp]
         case ("resp")
            charges = [ &
                      0.540724_wp, -0.597125_wp, 0.021522_wp, 0.017749_wp, 0.018962_wp, &
                      0.033311_wp, 0.019061_wp, 0.033003_wp, -0.25473_wp, 0.062684_wp, &
                      0.065789_wp, 0.081051_wp, -0.240021_wp, 0.058503_wp, 0.062135_wp, &
                      0.077383_wp]
         case ("mbis")
            charges = [ &
                      0.649284_wp, -0.571905_wp, -0.383145_wp, -0.383149_wp, 0.163714_wp, &
                      0.165088_wp, 0.163682_wp, 0.165115_wp, -0.392919_wp, 0.134404_wp, &
                      0.13476_wp, 0.139411_wp, -0.392912_wp, 0.134399_wp, 0.134765_wp, &
                      0.139413_wp]
         case ("chelpg")
            charges = [ &
                      0.541106_wp, -0.597423_wp, 0.02508_wp, 0.021024_wp, 0.018407_wp, &
                      0.032849_wp, 0.018566_wp, 0.032596_wp, -0.265391_wp, 0.065338_wp, &
                      0.068434_wp, 0.083857_wp, -0.250405_wp, 0.061103_wp, 0.064723_wp, &
                      0.080134_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (55)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070047_wp, -0.033363_wp, -0.041427_wp, 0.179883_wp, -0.261552_wp, &
                      -0.041888_wp, -0.033117_wp, -0.071663_wp, 0.021225_wp, 0.026784_wp, &
                      0.021725_wp, 0.025369_wp, 0.026369_wp, 0.035391_wp, 0.032502_wp, &
                      0.029158_wp, 0.037923_wp, 0.02218_wp, 0.026105_wp, 0.018462_wp, &
                      0.026462_wp, 0.023525_wp]
         case ("resp")
            charges = [ &
                      -0.378457_wp, 0.295639_wp, -0.284013_wp, 0.593621_wp, -0.520732_wp, &
                      -0.269299_wp, 0.25261_wp, -0.349909_wp, 0.086351_wp, 0.085274_wp, &
                      0.085823_wp, -0.030237_wp, -0.026914_wp, 0.072252_wp, 0.06457_wp, &
                      0.059643_wp, 0.060232_wp, -0.028329_wp, -0.018259_wp, 0.072734_wp, &
                      0.080063_wp, 0.097337_wp]
         case ("mbis")
            charges = [ &
                      -0.427128_wp, -0.124449_wp, -0.421228_wp, 0.619883_wp, -0.502398_wp, &
                      -0.42547_wp, -0.123719_wp, -0.423112_wp, 0.126813_wp, 0.135092_wp, &
                      0.127855_wp, 0.111624_wp, 0.116388_wp, 0.150421_wp, 0.148757_wp, &
                      0.145238_wp, 0.152335_wp, 0.097832_wp, 0.118858_wp, 0.120007_wp, &
                      0.132519_wp, 0.143888_wp]
         case ("chelpg")
            charges = [ &
                      -0.395885_wp, 0.32006_wp, -0.32033_wp, 0.617264_wp, -0.524674_wp, &
                      -0.29979_wp, 0.273323_wp, -0.368958_wp, 0.090544_wp, 0.088199_wp, &
                      0.090234_wp, -0.034344_wp, -0.030943_wp, 0.079406_wp, 0.072772_wp, &
                      0.067169_wp, 0.065952_wp, -0.030919_wp, -0.021773_wp, 0.07699_wp, &
                      0.083885_wp, 0.101817_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070521_wp, -0.03577_wp, -0.038424_wp, 0.18883_wp, -0.30561_wp, &
                      -0.038737_wp, -0.034407_wp, -0.074302_wp, 0.024392_wp, 0.026463_wp, &
                      0.024582_wp, 0.022072_wp, 0.020793_wp, 0.045765_wp, 0.043201_wp, &
                      0.039711_wp, 0.047613_wp, 0.027122_wp, 0.018959_wp, 0.023407_wp, &
                      0.026132_wp, 0.018735_wp]
         case ("resp")
            charges = [ &
                      -0.381476_wp, 0.324354_wp, -0.300094_wp, 0.654053_wp, -0.610529_wp, &
                      -0.29019_wp, 0.271695_wp, -0.383614_wp, 0.090432_wp, 0.081995_wp, &
                      0.089052_wp, -0.04417_wp, -0.044125_wp, 0.088949_wp, 0.081731_wp, &
                      0.079467_wp, 0.075614_wp, -0.020411_wp, -0.034135_wp, 0.089703_wp, &
                      0.084747_wp, 0.096951_wp]
         case ("mbis")
            charges = [ &
                      -0.436634_wp, -0.116602_wp, -0.447753_wp, 0.670586_wp, -0.575326_wp, &
                      -0.449767_wp, -0.118273_wp, -0.43683_wp, 0.133664_wp, 0.135661_wp, &
                      0.134156_wp, 0.105832_wp, 0.106518_wp, 0.170828_wp, 0.169853_wp, &
                      0.165797_wp, 0.170905_wp, 0.106759_wp, 0.106074_wp, 0.131778_wp, &
                      0.1347_wp, 0.138081_wp]
         case ("chelpg")
            charges = [ &
                      -0.399101_wp, 0.349636_wp, -0.33754_wp, 0.678403_wp, -0.614546_wp, &
                      -0.321705_wp, 0.292908_wp, -0.403086_wp, 0.094656_wp, 0.084915_wp, &
                      0.0935_wp, -0.048468_wp, -0.048346_wp, 0.096318_wp, 0.090171_wp, &
                      0.08726_wp, 0.081554_wp, -0.023038_wp, -0.037734_wp, 0.094059_wp, &
                      0.08865_wp, 0.101533_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070501_wp, -0.036079_wp, -0.038026_wp, 0.189745_wp, -0.310763_wp, &
                      -0.038302_wp, -0.034542_wp, -0.074574_wp, 0.024871_wp, 0.026313_wp, &
                      0.024993_wp, 0.021534_wp, 0.019943_wp, 0.047059_wp, 0.044596_wp, &
                      0.041132_wp, 0.04876_wp, 0.027798_wp, 0.01789_wp, 0.024156_wp, &
                      0.025998_wp, 0.018005_wp]
         case ("resp")
            charges = [ &
                      -0.381016_wp, 0.32651_wp, -0.300683_wp, 0.659754_wp, -0.62076_wp, &
                      -0.29149_wp, 0.273347_wp, -0.387243_wp, 0.090993_wp, 0.081259_wp, &
                      0.089456_wp, -0.04581_wp, -0.046271_wp, 0.090849_wp, 0.083816_wp, &
                      0.081937_wp, 0.077235_wp, -0.01914_wp, -0.036372_wp, 0.092019_wp, &
                      0.085105_wp, 0.096504_wp]
         case ("mbis")
            charges = [ &
                      -0.437477_wp, -0.115504_wp, -0.450715_wp, 0.676061_wp, -0.583793_wp, &
                      -0.452367_wp, -0.11751_wp, -0.438179_wp, 0.134604_wp, 0.135456_wp, &
                      0.134973_wp, 0.104867_wp, 0.104984_wp, 0.173304_wp, 0.172523_wp, &
                      0.168464_wp, 0.173029_wp, 0.107962_wp, 0.104108_wp, 0.133428_wp, &
                      0.13474_wp, 0.137048_wp]
         case ("chelpg")
            charges = [ &
                      -0.398657_wp, 0.351878_wp, -0.33824_wp, 0.684174_wp, -0.624783_wp, &
                      -0.323102_wp, 0.294605_wp, -0.406754_wp, 0.095219_wp, 0.084177_wp, &
                      0.093907_wp, -0.050128_wp, -0.050512_wp, 0.098239_wp, 0.092279_wp, &
                      0.089756_wp, 0.083196_wp, -0.021769_wp, -0.039979_wp, 0.096385_wp, &
                      0.089016_wp, 0.101095_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (56)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072873_wp, 0.181284_wp, -0.268737_wp, -0.042536_wp, -0.000963_wp, &
                      -0.071419_wp, -0.073271_wp, 0.041433_wp, 0.036729_wp, 0.042358_wp, &
                      0.037063_wp, 0.029778_wp, 0.025132_wp, 0.020327_wp, 0.025783_wp, &
                      0.022618_wp, 0.023965_wp, 0.02553_wp, 0.017803_wp]
         case ("resp")
            charges = [ &
                      -0.54562_wp, 0.72287_wp, -0.541436_wp, -0.502552_wp, 0.64587_wp, &
                      -0.559943_wp, -0.565043_wp, 0.155658_wp, 0.136796_wp, 0.144254_wp, &
                      0.113521_wp, 0.111971_wp, -0.060927_wp, 0.12645_wp, 0.124829_wp, &
                      0.114298_wp, 0.138121_wp, 0.124142_wp, 0.116742_wp]
         case ("mbis")
            charges = [ &
                      -0.638107_wp, 0.665333_wp, -0.517383_wp, -0.449062_wp, 0.078277_wp, &
                      -0.464599_wp, -0.459398_wp, 0.192479_wp, 0.174554_wp, 0.180343_wp, &
                      0.156147_wp, 0.151217_wp, 0.111122_wp, 0.13264_wp, 0.140715_wp, &
                      0.132535_wp, 0.150805_wp, 0.137331_wp, 0.125055_wp]
         case ("chelpg")
            charges = [ &
                      -0.572787_wp, 0.744795_wp, -0.545207_wp, -0.540561_wp, 0.67261_wp, &
                      -0.57987_wp, -0.587275_wp, 0.161663_wp, 0.143993_wp, 0.150469_wp, &
                      0.121318_wp, 0.120764_wp, -0.063513_wp, 0.131756_wp, 0.128656_wp, &
                      0.118974_wp, 0.143462_wp, 0.128781_wp, 0.121969_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071064_wp, 0.18981_wp, -0.315664_wp, -0.039462_wp, -0.001652_wp, &
                      -0.072095_wp, -0.076241_wp, 0.039235_wp, 0.0477_wp, 0.052149_wp, &
                      0.047425_wp, 0.040704_wp, 0.018457_wp, 0.023438_wp, 0.025742_wp, &
                      0.025302_wp, 0.018337_wp, 0.025282_wp, 0.0226_wp]
         case ("resp")
            charges = [ &
                      -0.568131_wp, 0.772426_wp, -0.629398_wp, -0.521123_wp, 0.697675_wp, &
                      -0.589132_wp, -0.616613_wp, 0.15427_wp, 0.159083_wp, 0.163733_wp, &
                      0.130323_wp, 0.132111_wp, -0.080983_wp, 0.135547_wp, 0.129254_wp, &
                      0.122009_wp, 0.139619_wp, 0.132461_wp, 0.136872_wp]
         case ("mbis")
            charges = [ &
                      -0.666122_wp, 0.717334_wp, -0.593745_wp, -0.476754_wp, 0.09528_wp, &
                      -0.478162_wp, -0.476771_wp, 0.191645_wp, 0.196356_wp, 0.19989_wp, &
                      0.176532_wp, 0.172965_wp, 0.098188_wp, 0.139926_wp, 0.142412_wp, &
                      0.139016_wp, 0.144274_wp, 0.140403_wp, 0.137339_wp]
         case ("chelpg")
            charges = [ &
                      -0.595589_wp, 0.79458_wp, -0.633197_wp, -0.55961_wp, 0.7248_wp, &
                      -0.609245_wp, -0.639093_wp, 0.160336_wp, 0.166362_wp, 0.170018_wp, &
                      0.13822_wp, 0.14101_wp, -0.083634_wp, 0.140896_wp, 0.133109_wp, &
                      0.126725_wp, 0.145016_wp, 0.137145_wp, 0.142152_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070819_wp, 0.190617_wp, -0.320784_wp, -0.039074_wp, -0.001756_wp, &
                      -0.072096_wp, -0.07655_wp, 0.038821_wp, 0.049023_wp, 0.053249_wp, &
                      0.048602_wp, 0.042079_wp, 0.017498_wp, 0.023928_wp, 0.025638_wp, &
                      0.025624_wp, 0.017563_wp, 0.025152_wp, 0.023289_wp]
         case ("resp")
            charges = [ &
                      -0.56979_wp, 0.77659_wp, -0.638759_wp, -0.521628_wp, 0.701894_wp, &
                      -0.591313_wp, -0.621799_wp, 0.153666_wp, 0.161623_wp, 0.165801_wp, &
                      0.13198_wp, 0.134423_wp, -0.083402_wp, 0.136667_wp, 0.129377_wp, &
                      0.122717_wp, 0.139418_wp, 0.133174_wp, 0.13936_wp]
         case ("mbis")
            charges = [ &
                      -0.66878_wp, 0.722576_wp, -0.602019_wp, -0.47959_wp, 0.097122_wp, &
                      -0.479348_wp, -0.478478_wp, 0.191176_wp, 0.198883_wp, 0.20201_wp, &
                      0.178774_wp, 0.175601_wp, 0.096351_wp, 0.140931_wp, 0.142345_wp, &
                      0.13971_wp, 0.143245_wp, 0.140529_wp, 0.138965_wp]
         case ("chelpg")
            charges = [ &
                      -0.597272_wp, 0.798764_wp, -0.642559_wp, -0.560158_wp, 0.729056_wp, &
                      -0.611442_wp, -0.6443_wp, 0.159737_wp, 0.168909_wp, 0.172092_wp, &
                      0.139885_wp, 0.143332_wp, -0.086059_wp, 0.142018_wp, 0.133234_wp, &
                      0.127436_wp, 0.14482_wp, 0.137862_wp, 0.144645_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (57)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.198133_wp, 0.035724_wp, -0.047728_wp, 0.027274_wp, -0.064215_wp, &
                      -0.047163_wp, 0.034547_wp, 0.039669_wp, 0.036271_wp, 0.036326_wp, &
                      0.03639_wp, 0.034093_wp, 0.037386_wp, 0.039555_wp]
         case ("resp")
            charges = [ &
                      -0.673639_wp, 0.486509_wp, -0.701807_wp, 0.742422_wp, -0.589056_wp, &
                      -0.687663_wp, 0.471015_wp, 0.035328_wp, 0.212885_wp, 0.150405_wp, &
                      0.150534_wp, 0.152208_wp, 0.21277_wp, 0.038089_wp]
         case ("mbis")
            charges = [ &
                      -0.45405_wp, 0.191134_wp, -0.330367_wp, 0.247666_wp, -0.479404_wp, &
                      -0.331795_wp, 0.186904_wp, 0.094521_wp, 0.157345_wp, 0.154737_wp, &
                      0.154819_wp, 0.151574_wp, 0.161614_wp, 0.095299_wp]
         case ("chelpg")
            charges = [ &
                      -0.685864_wp, 0.504874_wp, -0.728909_wp, 0.777035_wp, -0.61407_wp, &
                      -0.718342_wp, 0.491165_wp, 0.032548_wp, 0.218307_wp, 0.155548_wp, &
                      0.155642_wp, 0.157488_wp, 0.219639_wp, 0.034939_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.237901_wp, 0.02671_wp, -0.046509_wp, 0.031795_wp, -0.062775_wp, &
                      -0.045008_wp, 0.025274_wp, 0.044464_wp, 0.047524_wp, 0.041938_wp, &
                      0.041999_wp, 0.039565_wp, 0.048698_wp, 0.044222_wp]
         case ("resp")
            charges = [ &
                      -0.787367_wp, 0.521539_wp, -0.733797_wp, 0.776198_wp, -0.622907_wp, &
                      -0.717272_wp, 0.504651_wp, 0.040911_wp, 0.238307_wp, 0.165574_wp, &
                      0.16574_wp, 0.166625_wp, 0.238163_wp, 0.043635_wp]
         case ("mbis")
            charges = [ &
                      -0.530863_wp, 0.190953_wp, -0.341001_wp, 0.261816_wp, -0.498576_wp, &
                      -0.338932_wp, 0.185038_wp, 0.107381_wp, 0.178837_wp, 0.166016_wp, &
                      0.166089_wp, 0.162766_wp, 0.182315_wp, 0.108156_wp]
         case ("chelpg")
            charges = [ &
                      -0.799594_wp, 0.539889_wp, -0.760914_wp, 0.810923_wp, -0.648143_wp, &
                      -0.748007_wp, 0.524809_wp, 0.038141_wp, 0.243729_wp, 0.17077_wp, &
                      0.1709_wp, 0.171959_wp, 0.245048_wp, 0.040489_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.24273_wp, 0.02551_wp, -0.046158_wp, 0.032522_wp, -0.062493_wp, &
                      -0.044541_wp, 0.024048_wp, 0.044807_wp, 0.04888_wp, 0.042622_wp, &
                      0.042682_wp, 0.040219_wp, 0.050072_wp, 0.044553_wp]
         case ("resp")
            charges = [ &
                      -0.800984_wp, 0.525333_wp, -0.736582_wp, 0.779752_wp, -0.626141_wp, &
                      -0.719789_wp, 0.508298_wp, 0.041244_wp, 0.241109_wp, 0.167236_wp, &
                      0.167406_wp, 0.168174_wp, 0.240984_wp, 0.043962_wp]
         case ("mbis")
            charges = [ &
                      -0.540216_wp, 0.190778_wp, -0.341601_wp, 0.263497_wp, -0.50048_wp, &
                      -0.339137_wp, 0.184679_wp, 0.108554_wp, 0.181257_wp, 0.167294_wp, &
                      0.167366_wp, 0.16401_wp, 0.184668_wp, 0.109329_wp]
         case ("chelpg")
            charges = [ &
                      -0.813211_wp, 0.543684_wp, -0.763702_wp, 0.814489_wp, -0.651398_wp, &
                      -0.750531_wp, 0.528458_wp, 0.038473_wp, 0.246532_wp, 0.172436_wp, &
                      0.172571_wp, 0.173513_wp, 0.24787_wp, 0.040816_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (58)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.263976_wp, 0.178824_wp, -0.04133_wp, -0.035583_wp, -0.034216_wp, &
                      -0.069058_wp, -0.041326_wp, -0.035583_wp, -0.034218_wp, -0.069058_wp, &
                      0.033013_wp, 0.032909_wp, 0.02445_wp, 0.02448_wp, 0.018365_wp, &
                      0.018348_wp, 0.023843_wp, 0.023671_wp, 0.023679_wp, 0.03296_wp, &
                      0.03297_wp, 0.02447_wp, 0.024468_wp, 0.018354_wp, 0.018361_wp, &
                      0.023679_wp, 0.023673_wp, 0.023838_wp]
         case ("resp")
            charges = [ &
                      -0.583854_wp, 0.853176_wp, -0.522781_wp, 0.14913_wp, 0.149332_wp, &
                      -0.405221_wp, -0.518563_wp, 0.141663_wp, 0.154382_wp, -0.403778_wp, &
                      0.110937_wp, 0.110818_wp, 0.004914_wp, 0.004523_wp, -0.010618_wp, &
                      -0.010782_wp, 0.100028_wp, 0.092322_wp, 0.092241_wp, 0.11057_wp, &
                      0.11014_wp, 0.006816_wp, 0.006294_wp, -0.01239_wp, -0.012012_wp, &
                      0.090869_wp, 0.09192_wp, 0.099927_wp]
         case ("mbis")
            charges = [ &
                      -0.497602_wp, 0.611908_wp, -0.414441_wp, -0.177243_wp, -0.167542_wp, &
                      -0.421115_wp, -0.414436_wp, -0.177253_wp, -0.167553_wp, -0.42111_wp, &
                      0.146657_wp, 0.146532_wp, 0.119972_wp, 0.120116_wp, 0.099393_wp, &
                      0.099355_wp, 0.130081_wp, 0.130534_wp, 0.130543_wp, 0.146572_wp, &
                      0.146629_wp, 0.120065_wp, 0.120037_wp, 0.099381_wp, 0.099377_wp, &
                      0.130547_wp, 0.130532_wp, 0.130069_wp]
         case ("chelpg")
            charges = [ &
                      -0.588126_wp, 0.88313_wp, -0.562217_wp, 0.163763_wp, 0.159943_wp, &
                      -0.426133_wp, -0.55778_wp, 0.155891_wp, 0.165322_wp, -0.424897_wp, &
                      0.119341_wp, 0.119215_wp, 0.002229_wp, 0.001805_wp, -0.012337_wp, &
                      -0.012486_wp, 0.105242_wp, 0.096831_wp, 0.09678_wp, 0.118938_wp, &
                      0.118524_wp, 0.00419_wp, 0.003682_wp, -0.014153_wp, -0.013798_wp, &
                      0.095457_wp, 0.096475_wp, 0.105168_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.309758_wp, 0.188169_wp, -0.038428_wp, -0.03807_wp, -0.033016_wp, &
                      -0.071739_wp, -0.038425_wp, -0.038071_wp, -0.033017_wp, -0.071739_wp, &
                      0.043647_wp, 0.043552_wp, 0.02002_wp, 0.019971_wp, 0.021814_wp, &
                      0.021806_wp, 0.024673_wp, 0.023278_wp, 0.023277_wp, 0.043602_wp, &
                      0.043607_wp, 0.019993_wp, 0.020009_wp, 0.02181_wp, 0.021813_wp, &
                      0.02328_wp, 0.023278_wp, 0.024668_wp]
         case ("resp")
            charges = [ &
                      -0.678644_wp, 0.92973_wp, -0.553593_wp, 0.163975_wp, 0.185122_wp, &
                      -0.442373_wp, -0.548789_wp, 0.155675_wp, 0.190748_wp, -0.440939_wp, &
                      0.130244_wp, 0.130094_wp, -0.007869_wp, -0.008512_wp, -0.013146_wp, &
                      -0.013282_wp, 0.107672_wp, 0.097811_wp, 0.097743_wp, 0.129757_wp, &
                      0.129359_wp, -0.005951_wp, -0.00641_wp, -0.015024_wp, -0.014677_wp, &
                      0.096329_wp, 0.097415_wp, 0.107534_wp]
         case ("mbis")
            charges = [ &
                      -0.573042_wp, 0.665226_wp, -0.442419_wp, -0.172157_wp, -0.16555_wp, &
                      -0.436057_wp, -0.442415_wp, -0.172174_wp, -0.165557_wp, -0.436053_wp, &
                      0.167838_wp, 0.167731_wp, 0.112474_wp, 0.112477_wp, 0.104903_wp, &
                      0.104878_wp, 0.134414_wp, 0.13269_wp, 0.132682_wp, 0.167768_wp, &
                      0.167815_wp, 0.112484_wp, 0.112487_wp, 0.1049_wp, 0.10489_wp, &
                      0.132691_wp, 0.132682_wp, 0.134401_wp]
         case ("chelpg")
            charges = [ &
                      -0.682958_wp, 0.960267_wp, -0.593937_wp, 0.179312_wp, 0.196113_wp, &
                      -0.463721_wp, -0.588927_wp, 0.17062_wp, 0.202053_wp, -0.46249_wp, &
                      0.138825_wp, 0.138668_wp, -0.010749_wp, -0.011426_wp, -0.014983_wp, &
                      -0.015104_wp, 0.112992_wp, 0.1024_wp, 0.102363_wp, 0.138304_wp, &
                      0.137922_wp, -0.008772_wp, -0.009218_wp, -0.016901_wp, -0.016578_wp, &
                      0.100996_wp, 0.102049_wp, 0.112881_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.31602_wp, 0.189293_wp, -0.037963_wp, -0.03843_wp, -0.032776_wp, &
                      -0.072076_wp, -0.037959_wp, -0.038431_wp, -0.032777_wp, -0.072076_wp, &
                      0.045246_wp, 0.045154_wp, 0.019198_wp, 0.019136_wp, 0.022417_wp, &
                      0.022411_wp, 0.02483_wp, 0.023105_wp, 0.023102_wp, 0.045202_wp, &
                      0.045207_wp, 0.019163_wp, 0.019181_wp, 0.022415_wp, 0.022417_wp, &
                      0.023105_wp, 0.023104_wp, 0.024824_wp]
         case ("resp")
            charges = [ &
                      -0.691096_wp, 0.937688_wp, -0.555412_wp, 0.164587_wp, 0.190397_wp, &
                      -0.447225_wp, -0.550505_wp, 0.156171_wp, 0.196101_wp, -0.445789_wp, &
                      0.132743_wp, 0.132591_wp, -0.009756_wp, -0.010438_wp, -0.013296_wp, &
                      -0.013426_wp, 0.108773_wp, 0.098325_wp, 0.098258_wp, 0.132232_wp, &
                      0.131842_wp, -0.007839_wp, -0.008289_wp, -0.015187_wp, -0.014846_wp, &
                      0.096837_wp, 0.097929_wp, 0.108629_wp]
         case ("mbis")
            charges = [ &
                      -0.58329_wp, 0.671984_wp, -0.446062_wp, -0.171112_wp, -0.165374_wp, &
                      -0.437752_wp, -0.446059_wp, -0.17113_wp, -0.165381_wp, -0.437749_wp, &
                      0.170921_wp, 0.17082_wp, 0.111015_wp, 0.110996_wp, 0.105875_wp, &
                      0.105855_wp, 0.135038_wp, 0.13272_wp, 0.132708_wp, 0.170855_wp, &
                      0.170901_wp, 0.111013_wp, 0.11102_wp, 0.105876_wp, 0.105863_wp, &
                      0.132719_wp, 0.132712_wp, 0.135025_wp]
         case ("chelpg")
            charges = [ &
                      -0.695411_wp, 0.968277_wp, -0.595844_wp, 0.180001_wp, 0.201434_wp, &
                      -0.468622_wp, -0.590731_wp, 0.171192_wp, 0.207452_wp, -0.467389_wp, &
                      0.141341_wp, 0.141183_wp, -0.012659_wp, -0.013376_wp, -0.015147_wp, &
                      -0.015263_wp, 0.114105_wp, 0.102922_wp, 0.102886_wp, 0.140797_wp, &
                      0.140423_wp, -0.010684_wp, -0.011121_wp, -0.017078_wp, -0.016761_wp, &
                      0.101513_wp, 0.102572_wp, 0.113987_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (59)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.234826_wp, -0.306326_wp, -0.190053_wp, 0.179894_wp, -0.062249_wp, &
                      0.04752_wp, 0.048629_wp, 0.047757_wp]
         case ("resp")
            charges = [ &
                      0.741826_wp, -0.554496_wp, -0.57546_wp, 0.406948_wp, -0.474671_wp, &
                      0.144435_wp, 0.155519_wp, 0.1559_wp]
         case ("mbis")
            charges = [ &
                      0.871815_wp, -0.613281_wp, -0.653034_wp, 0.460601_wp, -0.649238_wp, &
                      0.198258_wp, 0.192478_wp, 0.192401_wp]
         case ("chelpg")
            charges = [ &
                      0.748727_wp, -0.556023_wp, -0.576399_wp, 0.40674_wp, -0.49114_wp, &
                      0.14868_wp, 0.159499_wp, 0.159915_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.240183_wp, -0.342241_wp, -0.190583_wp, 0.194431_wp, -0.061907_wp, &
                      0.049106_wp, 0.05607_wp, 0.054939_wp]
         case ("resp")
            charges = [ &
                      0.789275_wp, -0.619759_wp, -0.607579_wp, 0.445092_wp, -0.507572_wp, &
                      0.153111_wp, 0.173636_wp, 0.173796_wp]
         case ("mbis")
            charges = [ &
                      0.913488_wp, -0.670797_wp, -0.6774_wp, 0.491538_wp, -0.677899_wp, &
                      0.204404_wp, 0.208584_wp, 0.208081_wp]
         case ("chelpg")
            charges = [ &
                      0.796236_wp, -0.621297_wp, -0.608513_wp, 0.444871_wp, -0.524183_wp, &
                      0.157395_wp, 0.177647_wp, 0.177844_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.241252_wp, -0.35046_wp, -0.1902_wp, 0.197434_wp, -0.061758_wp, &
                      0.04924_wp, 0.057834_wp, 0.056655_wp]
         case ("resp")
            charges = [ &
                      0.798785_wp, -0.634493_wp, -0.613241_wp, 0.452948_wp, -0.514308_wp, &
                      0.154528_wp, 0.17783_wp, 0.17795_wp]
         case ("mbis")
            charges = [ &
                      0.92216_wp, -0.683786_wp, -0.681442_wp, 0.497628_wp, -0.683864_wp, &
                      0.205327_wp, 0.212279_wp, 0.211696_wp]
         case ("chelpg")
            charges = [ &
                      0.805764_wp, -0.636043_wp, -0.614175_wp, 0.45273_wp, -0.530964_wp, &
                      0.158825_wp, 0.181854_wp, 0.18201_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (60)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.280608_wp, 0.180428_wp, -0.071524_wp, -0.071545_wp, 0.037426_wp, &
                      0.041975_wp, 0.042221_wp, 0.037559_wp, 0.042195_wp, 0.041875_wp]
         case ("resp")
            charges = [ &
                      -0.542928_wp, 0.772041_wp, -0.576499_wp, -0.575259_wp, 0.145591_wp, &
                      0.152954_wp, 0.162194_wp, 0.146155_wp, 0.162422_wp, 0.153328_wp]
         case ("mbis")
            charges = [ &
                      -0.523545_wp, 0.695011_wp, -0.628796_wp, -0.628924_wp, 0.174051_wp, &
                      0.177788_wp, 0.191332_wp, 0.174159_wp, 0.191293_wp, 0.177632_wp]
         case ("chelpg")
            charges = [ &
                      -0.545285_wp, 0.786523_wp, -0.597855_wp, -0.596742_wp, 0.15097_wp, &
                      0.15788_wp, 0.167298_wp, 0.151499_wp, 0.167487_wp, 0.158224_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.336063_wp, 0.189442_wp, -0.069325_wp, -0.069338_wp, 0.049694_wp, &
                      0.053396_wp, 0.039541_wp, 0.049832_wp, 0.039487_wp, 0.053334_wp]
         case ("resp")
            charges = [ &
                      -0.634909_wp, 0.814475_wp, -0.591664_wp, -0.590332_wp, 0.168963_wp, &
                      0.172995_wp, 0.158627_wp, 0.169538_wp, 0.15886_wp, 0.173446_wp]
         case ("mbis")
            charges = [ &
                      -0.610666_wp, 0.751676_wp, -0.65812_wp, -0.658185_wp, 0.197963_wp, &
                      0.200101_wp, 0.1896_wp, 0.198083_wp, 0.189524_wp, 0.200025_wp]
         case ("chelpg")
            charges = [ &
                      -0.637278_wp, 0.829062_wp, -0.613194_wp, -0.611991_wp, 0.174387_wp, &
                      0.177963_wp, 0.163773_wp, 0.174927_wp, 0.163968_wp, 0.178384_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.33968_wp, 0.189948_wp, -0.069152_wp, -0.069164_wp, 0.050577_wp, &
                      0.054182_wp, 0.039255_wp, 0.050714_wp, 0.039199_wp, 0.054123_wp]
         case ("resp")
            charges = [ &
                      -0.64077_wp, 0.816432_wp, -0.592054_wp, -0.590714_wp, 0.170564_wp, &
                      0.174287_wp, 0.158068_wp, 0.17114_wp, 0.158303_wp, 0.174744_wp]
         case ("mbis")
            charges = [ &
                      -0.616312_wp, 0.755086_wp, -0.6598_wp, -0.65986_wp, 0.199623_wp, &
                      0.201586_wp, 0.189248_wp, 0.199743_wp, 0.18917_wp, 0.201517_wp]
         case ("chelpg")
            charges = [ &
                      -0.643139_wp, 0.831019_wp, -0.613589_wp, -0.612377_wp, 0.175991_wp, &
                      0.179256_wp, 0.163214_wp, 0.176531_wp, 0.16341_wp, 0.179683_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (61)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.24597_wp, 0.084299_wp, -0.024532_wp, 0.062071_wp, 0.062067_wp, &
                      0.062064_wp]
         case ("resp")
            charges = [ &
                      -0.489044_wp, 0.448513_wp, -0.525349_wp, 0.188778_wp, 0.188306_wp, &
                      0.188796_wp]
         case ("mbis")
            charges = [ &
                      -0.467955_wp, 0.400633_wp, -0.507861_wp, 0.191724_wp, 0.191738_wp, &
                      0.191721_wp]
         case ("chelpg")
            charges = [ &
                      -0.490949_wp, 0.455488_wp, -0.544663_wp, 0.193521_wp, 0.193063_wp, &
                      0.193539_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.320352_wp, 0.107261_wp, -0.015943_wp, 0.07636_wp, 0.0763_wp, &
                      0.076374_wp]
         case ("resp")
            charges = [ &
                      -0.607973_wp, 0.51968_wp, -0.555394_wp, 0.214689_wp, 0.214225_wp, &
                      0.214773_wp]
         case ("mbis")
            charges = [ &
                      -0.597899_wp, 0.496606_wp, -0.554302_wp, 0.21855_wp, 0.218463_wp, &
                      0.218582_wp]
         case ("chelpg")
            charges = [ &
                      -0.609897_wp, 0.526727_wp, -0.574918_wp, 0.219484_wp, 0.219034_wp, &
                      0.21957_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.323126_wp, 0.108108_wp, -0.015606_wp, 0.07689_wp, 0.076828_wp, &
                      0.076905_wp]
         case ("resp")
            charges = [ &
                      -0.612419_wp, 0.52227_wp, -0.556311_wp, 0.215612_wp, 0.215148_wp, &
                      0.2157_wp]
         case ("mbis")
            charges = [ &
                      -0.602728_wp, 0.500142_wp, -0.555936_wp, 0.219527_wp, 0.219435_wp, &
                      0.21956_wp]
         case ("chelpg")
            charges = [ &
                      -0.614342_wp, 0.529316_wp, -0.575836_wp, 0.220408_wp, 0.219958_wp, &
                      0.220496_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (62)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.022221_wp, -0.035135_wp, -0.028032_wp, -0.017679_wp, -0.016944_wp, &
                      -0.030628_wp, 0.040823_wp, 0.038537_wp, 0.036802_wp, 0.171172_wp, &
                      0.042331_wp, 0.040481_wp, -0.270522_wp, -0.071224_wp, 0.042368_wp, &
                      0.039957_wp, 0.039911_wp]
         case ("resp")
            charges = [ &
                      -0.10462_wp, -0.134271_wp, -0.151107_wp, -0.039142_wp, -0.116894_wp, &
                      -0.125856_wp, 0.127959_wp, 0.129809_wp, 0.128903_wp, 0.593386_wp, &
                      0.137503_wp, 0.128946_wp, -0.504348_wp, -0.471401_wp, 0.145179_wp, &
                      0.127959_wp, 0.127991_wp]
         case ("mbis")
            charges = [ &
                      -0.099066_wp, -0.157214_wp, -0.113458_wp, -0.133147_wp, -0.083226_wp, &
                      -0.156364_wp, 0.135814_wp, 0.141029_wp, 0.140186_wp, 0.613344_wp, &
                      0.153056_wp, 0.142764_wp, -0.508103_wp, -0.618361_wp, 0.192259_wp, &
                      0.175281_wp, 0.175203_wp]
         case ("chelpg")
            charges = [ &
                      -0.106533_wp, -0.131713_wp, -0.152902_wp, -0.041509_wp, -0.119124_wp, &
                      -0.12322_wp, 0.128212_wp, 0.12933_wp, 0.129463_wp, 0.604578_wp, &
                      0.138048_wp, 0.128458_wp, -0.506664_wp, -0.493623_wp, 0.150698_wp, &
                      0.133236_wp, 0.133264_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.022458_wp, -0.035885_wp, -0.02447_wp, -0.024297_wp, -0.027079_wp, &
                      -0.035113_wp, 0.050567_wp, 0.048359_wp, 0.049725_wp, 0.176183_wp, &
                      0.038615_wp, 0.048215_wp, -0.317638_wp, -0.069155_wp, 0.040643_wp, &
                      0.051908_wp, 0.051874_wp]
         case ("resp")
            charges = [ &
                      -0.114929_wp, -0.140218_wp, -0.167833_wp, -0.040241_wp, -0.132806_wp, &
                      -0.145475_wp, 0.148522_wp, 0.148572_wp, 0.157498_wp, 0.620981_wp, &
                      0.140803_wp, 0.149471_wp, -0.58347_wp, -0.484714_wp, 0.1432_wp, &
                      0.150323_wp, 0.150317_wp]
         case ("mbis")
            charges = [ &
                      -0.105078_wp, -0.17267_wp, -0.109701_wp, -0.15447_wp, -0.098433_wp, &
                      -0.171047_wp, 0.154716_wp, 0.160929_wp, 0.16277_wp, 0.658651_wp, &
                      0.153366_wp, 0.160348_wp, -0.581672_wp, -0.648115_wp, 0.192502_wp, &
                      0.198976_wp, 0.198922_wp]
         case ("chelpg")
            charges = [ &
                      -0.117068_wp, -0.137389_wp, -0.169953_wp, -0.042301_wp, -0.135286_wp, &
                      -0.142648_wp, 0.148819_wp, 0.148039_wp, 0.15813_wp, 0.632164_wp, &
                      0.141395_wp, 0.14896_wp, -0.585781_wp, -0.507153_wp, 0.14878_wp, &
                      0.155652_wp, 0.155641_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.022308_wp, -0.0357_wp, -0.023995_wp, -0.024862_wp, -0.028041_wp, &
                      -0.035476_wp, 0.051329_wp, 0.049179_wp, 0.050813_wp, 0.17644_wp, &
                      0.038081_wp, 0.048714_wp, -0.32137_wp, -0.068973_wp, 0.040392_wp, &
                      0.052905_wp, 0.052872_wp]
         case ("resp")
            charges = [ &
                      -0.115372_wp, -0.140128_wp, -0.168947_wp, -0.040187_wp, -0.134365_wp, &
                      -0.146942_wp, 0.150002_wp, 0.14994_wp, 0.159822_wp, 0.622429_wp, &
                      0.140709_wp, 0.150884_wp, -0.589564_wp, -0.48528_wp, 0.142758_wp, &
                      0.152125_wp, 0.152115_wp]
         case ("mbis")
            charges = [ &
                      -0.105213_wp, -0.173453_wp, -0.109029_wp, -0.156252_wp, -0.099828_wp, &
                      -0.172122_wp, 0.156084_wp, 0.162426_wp, 0.164553_wp, 0.661903_wp, &
                      0.153031_wp, 0.161525_wp, -0.587426_wp, -0.650238_wp, 0.192283_wp, &
                      0.200902_wp, 0.20085_wp]
         case ("chelpg")
            charges = [ &
                      -0.117539_wp, -0.137264_wp, -0.171103_wp, -0.042218_wp, -0.136883_wp, &
                      -0.144089_wp, 0.150306_wp, 0.149404_wp, 0.160463_wp, 0.633617_wp, &
                      0.141309_wp, 0.150372_wp, -0.591887_wp, -0.50773_wp, 0.148342_wp, &
                      0.157457_wp, 0.157443_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (63)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.059934_wp, -0.034549_wp, -0.062496_wp, 0.057555_wp, -0.062492_wp, &
                      -0.034548_wp, 0.030367_wp, 0.034979_wp, 0.029725_wp, -0.191412_wp, &
                      0.029728_wp, 0.034978_wp, 0.11405_wp, 0.114046_wp]
         case ("resp")
            charges = [ &
                      -0.165072_wp, -0.118661_wp, -0.283282_wp, 0.390096_wp, -0.27913_wp, &
                      -0.123229_wp, 0.123456_wp, 0.131293_wp, 0.157045_wp, -0.767785_wp, &
                      0.156437_wp, 0.132611_wp, 0.322616_wp, 0.323603_wp]
         case ("mbis")
            charges = [ &
                      -0.222922_wp, -0.069359_wp, -0.317153_wp, 0.416365_wp, -0.317135_wp, &
                      -0.069362_wp, 0.141652_wp, 0.128137_wp, 0.147301_wp, -0.854921_wp, &
                      0.147305_wp, 0.128135_wp, 0.370978_wp, 0.370976_wp]
         case ("chelpg")
            charges = [ &
                      -0.178446_wp, -0.104511_wp, -0.300755_wp, 0.409289_wp, -0.296577_wp, &
                      -0.109137_wp, 0.125796_wp, 0.12907_wp, 0.160442_wp, -0.772594_wp, &
                      0.159827_wp, 0.130413_wp, 0.323101_wp, 0.324085_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071494_wp, -0.043899_wp, -0.069227_wp, 0.051183_wp, -0.069225_wp, &
                      -0.043899_wp, 0.033846_wp, 0.040212_wp, 0.038371_wp, -0.197826_wp, &
                      0.038373_wp, 0.040211_wp, 0.126688_wp, 0.126683_wp]
         case ("resp")
            charges = [ &
                      -0.211177_wp, -0.125077_wp, -0.337765_wp, 0.474823_wp, -0.333412_wp, &
                      -0.129828_wp, 0.143089_wp, 0.147705_wp, 0.182599_wp, -0.908044_wp, &
                      0.181867_wp, 0.1491_wp, 0.382637_wp, 0.383483_wp]
         case ("mbis")
            charges = [ &
                      -0.24871_wp, -0.087427_wp, -0.341339_wp, 0.426967_wp, -0.341321_wp, &
                      -0.08743_wp, 0.155675_wp, 0.143654_wp, 0.168978_wp, -0.906347_wp, &
                      0.16898_wp, 0.143651_wp, 0.402335_wp, 0.402331_wp]
         case ("chelpg")
            charges = [ &
                      -0.225033_wp, -0.110482_wp, -0.35568_wp, 0.494457_wp, -0.351305_wp, &
                      -0.115289_wp, 0.145529_wp, 0.145406_wp, 0.186076_wp, -0.912939_wp, &
                      0.185338_wp, 0.146825_wp, 0.383128_wp, 0.38397_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074043_wp, -0.045881_wp, -0.070467_wp, 0.05004_wp, -0.070465_wp, &
                      -0.045881_wp, 0.034346_wp, 0.041148_wp, 0.040169_wp, -0.198833_wp, &
                      0.04017_wp, 0.041146_wp, 0.129277_wp, 0.129272_wp]
         case ("resp")
            charges = [ &
                      -0.221729_wp, -0.126006_wp, -0.349159_wp, 0.494268_wp, -0.344752_wp, &
                      -0.130809_wp, 0.147049_wp, 0.150857_wp, 0.187744_wp, -0.939172_wp, &
                      0.186981_wp, 0.152273_wp, 0.395821_wp, 0.396633_wp]
         case ("mbis")
            charges = [ &
                      -0.254234_wp, -0.091187_wp, -0.345804_wp, 0.428952_wp, -0.345787_wp, &
                      -0.09119_wp, 0.158258_wp, 0.146609_wp, 0.17333_wp, -0.915675_wp, &
                      0.173331_wp, 0.146606_wp, 0.408396_wp, 0.408392_wp]
         case ("chelpg")
            charges = [ &
                      -0.235682_wp, -0.111324_wp, -0.367167_wp, 0.514_wp, -0.362737_wp, &
                      -0.116182_wp, 0.149507_wp, 0.148545_wp, 0.191237_wp, -0.944092_wp, &
                      0.190467_wp, 0.149985_wp, 0.396318_wp, 0.397126_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (64)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.053397_wp, -0.031758_wp, -0.052608_wp, 0.082872_wp, -0.064592_wp, &
                      -0.033996_wp, 0.032725_wp, 0.037257_wp, -0.157558_wp, 0.031198_wp, &
                      0.035839_wp, 0.026173_wp, 0.039216_wp, 0.032453_wp, 0.032463_wp, &
                      0.043714_wp]
         case ("resp")
            charges = [ &
                      -0.23886_wp, -0.01394_wp, -0.390682_wp, 0.485964_wp, -0.273818_wp, &
                      -0.117594_wp, 0.144218_wp, 0.119382_wp, -0.336829_wp, 0.137261_wp, &
                      0.14015_wp, -0.108829_wp, 0.176156_wp, 0.079191_wp, 0.079835_wp, &
                      0.118394_wp]
         case ("mbis")
            charges = [ &
                      -0.210055_wp, -0.072343_wp, -0.301612_wp, 0.391895_wp, -0.351777_wp, &
                      -0.061848_wp, 0.14225_wp, 0.130694_wp, -0.327822_wp, 0.163394_wp, &
                      0.128221_wp, -0.109616_wp, 0.167315_wp, 0.095687_wp, 0.095733_wp, &
                      0.119883_wp]
         case ("chelpg")
            charges = [ &
                      -0.259395_wp, 0.008242_wp, -0.416975_wp, 0.513993_wp, -0.298931_wp, &
                      -0.097289_wp, 0.147916_wp, 0.115594_wp, -0.340256_wp, 0.142486_wp, &
                      0.136956_wp, -0.119058_wp, 0.181324_wp, 0.08147_wp, 0.082116_wp, &
                      0.121807_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.059685_wp, -0.037507_wp, -0.06089_wp, 0.076235_wp, -0.066167_wp, &
                      -0.038543_wp, 0.037675_wp, 0.043012_wp, -0.168727_wp, 0.039995_wp, &
                      0.042174_wp, 0.026684_wp, 0.041586_wp, 0.038355_wp, 0.03843_wp, &
                      0.047374_wp]
         case ("resp")
            charges = [ &
                      -0.264042_wp, -0.023921_wp, -0.414892_wp, 0.498615_wp, -0.292809_wp, &
                      -0.125709_wp, 0.160339_wp, 0.134973_wp, -0.37407_wp, 0.159023_wp, &
                      0.154975_wp, -0.105897_wp, 0.18811_wp, 0.090118_wp, 0.090843_wp, &
                      0.124346_wp]
         case ("mbis")
            charges = [ &
                      -0.227785_wp, -0.082352_wp, -0.324271_wp, 0.39103_wp, -0.360904_wp, &
                      -0.075186_wp, 0.155953_wp, 0.144711_wp, -0.344912_wp, 0.181599_wp, &
                      0.143211_wp, -0.122391_wp, 0.177913_wp, 0.107675_wp, 0.10784_wp, &
                      0.127868_wp]
         case ("chelpg")
            charges = [ &
                      -0.283907_wp, -0.002525_wp, -0.440543_wp, 0.526202_wp, -0.317468_wp, &
                      -0.105922_wp, 0.163904_wp, 0.131363_wp, -0.377442_wp, 0.164172_wp, &
                      0.151869_wp, -0.116174_wp, 0.193148_wp, 0.092411_wp, 0.093139_wp, &
                      0.127772_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.061836_wp, -0.039733_wp, -0.064385_wp, 0.073574_wp, -0.066337_wp, &
                      -0.03977_wp, 0.039478_wp, 0.044963_wp, -0.173018_wp, 0.043555_wp, &
                      0.044672_wp, 0.026938_wp, 0.041969_wp, 0.040657_wp, 0.04076_wp, &
                      0.048513_wp]
         case ("resp")
            charges = [ &
                      -0.273067_wp, -0.027793_wp, -0.42415_wp, 0.502455_wp, -0.299236_wp, &
                      -0.127536_wp, 0.166126_wp, 0.140522_wp, -0.388389_wp, 0.167517_wp, &
                      0.160346_wp, -0.104069_wp, 0.191889_wp, 0.094304_wp, 0.095068_wp, &
                      0.126013_wp]
         case ("mbis")
            charges = [ &
                      -0.234077_wp, -0.085921_wp, -0.333271_wp, 0.390177_wp, -0.363407_wp, &
                      -0.079315_wp, 0.160852_wp, 0.149586_wp, -0.351475_wp, 0.188589_wp, &
                      0.148738_wp, -0.12674_wp, 0.181182_wp, 0.112248_wp, 0.112466_wp, &
                      0.130369_wp]
         case ("chelpg")
            charges = [ &
                      -0.292679_wp, -0.00669_wp, -0.449548_wp, 0.52988_wp, -0.323742_wp, &
                      -0.107954_wp, 0.169646_wp, 0.136978_wp, -0.391752_wp, 0.172645_wp, &
                      0.157276_wp, -0.114345_wp, 0.196882_wp, 0.096598_wp, 0.097365_wp, &
                      0.129439_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (65)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.017798_wp, -0.032127_wp, -0.021681_wp, -0.016363_wp, -0.012319_wp, &
                      -0.028249_wp, 0.042481_wp, 0.040849_wp, 0.040103_wp, 0.144704_wp, &
                      0.046251_wp, 0.042045_wp, -0.260543_wp, 0.032651_wp]
         case ("resp")
            charges = [ &
                      -0.135539_wp, -0.072146_wp, -0.239239_wp, 0.100257_wp, -0.116668_wp, &
                      -0.127682_wp, 0.136325_wp, 0.124298_wp, 0.144034_wp, 0.314856_wp, &
                      0.135755_wp, 0.134156_wp, -0.431721_wp, 0.033313_wp]
         case ("mbis")
            charges = [ &
                      -0.094632_wp, -0.149444_wp, -0.116587_wp, -0.102308_wp, -0.079181_wp, &
                      -0.157332_wp, 0.137313_wp, 0.143068_wp, 0.134258_wp, 0.415689_wp, &
                      0.152424_wp, 0.14506_wp, -0.460744_wp, 0.03242_wp]
         case ("chelpg")
            charges = [ &
                      -0.146302_wp, -0.060301_wp, -0.253372_wp, 0.114268_wp, -0.129764_wp, &
                      -0.11689_wp, 0.138214_wp, 0.122316_wp, 0.146861_wp, 0.312669_wp, &
                      0.138298_wp, 0.132507_wp, -0.431845_wp, 0.033341_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.016466_wp, -0.032879_wp, -0.018236_wp, -0.023542_wp, -0.020361_wp, &
                      -0.031387_wp, 0.052711_wp, 0.050477_wp, 0.053297_wp, 0.151662_wp, &
                      0.044774_wp, 0.050345_wp, -0.309914_wp, 0.049524_wp]
         case ("resp")
            charges = [ &
                      -0.142095_wp, -0.079177_wp, -0.257472_wp, 0.09592_wp, -0.136132_wp, &
                      -0.141156_wp, 0.156314_wp, 0.142597_wp, 0.172897_wp, 0.348695_wp, &
                      0.14297_wp, 0.153682_wp, -0.516264_wp, 0.059222_wp]
         case ("mbis")
            charges = [ &
                      -0.097553_wp, -0.164913_wp, -0.115365_wp, -0.123636_wp, -0.091692_wp, &
                      -0.170888_wp, 0.156227_wp, 0.162586_wp, 0.157317_wp, 0.444269_wp, &
                      0.155417_wp, 0.162943_wp, -0.534731_wp, 0.060023_wp]
         case ("chelpg")
            charges = [ &
                      -0.152758_wp, -0.067472_wp, -0.271552_wp, 0.110044_wp, -0.149266_wp, &
                      -0.130424_wp, 0.158185_wp, 0.140655_wp, 0.175704_wp, 0.346345_wp, &
                      0.145534_wp, 0.152048_wp, -0.516338_wp, 0.059297_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.016204_wp, -0.032734_wp, -0.017827_wp, -0.024113_wp, -0.021104_wp, &
                      -0.031603_wp, 0.053476_wp, 0.051221_wp, 0.054308_wp, 0.152078_wp, &
                      0.044442_wp, 0.050882_wp, -0.313582_wp, 0.050765_wp]
         case ("resp")
            charges = [ &
                      -0.142135_wp, -0.079431_wp, -0.258185_wp, 0.095024_wp, -0.137642_wp, &
                      -0.142065_wp, 0.157673_wp, 0.143875_wp, 0.174924_wp, 0.351328_wp, &
                      0.143127_wp, 0.15496_wp, -0.522556_wp, 0.061104_wp]
         case ("mbis")
            charges = [ &
                      -0.09744_wp, -0.16573_wp, -0.114919_wp, -0.125353_wp, -0.092793_wp, &
                      -0.171775_wp, 0.157537_wp, 0.163966_wp, 0.158983_wp, 0.446248_wp, &
                      0.155308_wp, 0.164103_wp, -0.540195_wp, 0.062066_wp]
         case ("chelpg")
            charges = [ &
                      -0.15278_wp, -0.067741_wp, -0.272251_wp, 0.109126_wp, -0.150767_wp, &
                      -0.131348_wp, 0.159542_wp, 0.141939_wp, 0.17773_wp, 0.348991_wp, &
                      0.145691_wp, 0.153332_wp, -0.522646_wp, 0.061181_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (66)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.035681_wp, -0.035676_wp, -0.035678_wp, -0.035682_wp, -0.035677_wp, &
                      -0.035678_wp, 0.03568_wp, 0.035678_wp, 0.03568_wp, 0.035678_wp, &
                      0.035677_wp, 0.035681_wp]
         case ("resp")
            charges = [ &
                      -0.11978_wp, -0.121266_wp, -0.117473_wp, -0.121759_wp, -0.120002_wp, &
                      -0.122077_wp, 0.120308_wp, 0.12062_wp, 0.11903_wp, 0.120734_wp, &
                      0.120588_wp, 0.121078_wp]
         case ("mbis")
            charges = [ &
                      -0.134219_wp, -0.134217_wp, -0.134202_wp, -0.134186_wp, -0.134211_wp, &
                      -0.134225_wp, 0.134219_wp, 0.134206_wp, 0.13421_wp, 0.134206_wp, &
                      0.134203_wp, 0.134217_wp]
         case ("chelpg")
            charges = [ &
                      -0.119757_wp, -0.121528_wp, -0.117427_wp, -0.122023_wp, -0.119981_wp, &
                      -0.122322_wp, 0.120397_wp, 0.120761_wp, 0.119111_wp, 0.120876_wp, &
                      0.120678_wp, 0.121215_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.039597_wp, -0.03959_wp, -0.039583_wp, -0.039597_wp, -0.03959_wp, &
                      -0.039584_wp, 0.039588_wp, 0.039586_wp, 0.039598_wp, 0.039586_wp, &
                      0.039586_wp, 0.039599_wp]
         case ("resp")
            charges = [ &
                      -0.130217_wp, -0.131755_wp, -0.127792_wp, -0.132262_wp, -0.130418_wp, &
                      -0.132492_wp, 0.130737_wp, 0.131052_wp, 0.129447_wp, 0.131169_wp, &
                      0.131013_wp, 0.131516_wp]
         case ("mbis")
            charges = [ &
                      -0.144148_wp, -0.144127_wp, -0.144089_wp, -0.14411_wp, -0.144121_wp, &
                      -0.144115_wp, 0.144126_wp, 0.144109_wp, 0.144126_wp, 0.144112_wp, &
                      0.144106_wp, 0.144133_wp]
         case ("chelpg")
            charges = [ &
                      -0.130214_wp, -0.132003_wp, -0.127767_wp, -0.132511_wp, -0.130418_wp, &
                      -0.132724_wp, 0.130833_wp, 0.131194_wp, 0.129536_wp, 0.131311_wp, &
                      0.13111_wp, 0.131653_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.043129_wp, -0.043118_wp, -0.043102_wp, -0.043127_wp, -0.043118_wp, &
                      -0.043104_wp, 0.043111_wp, 0.043109_wp, 0.043131_wp, 0.043108_wp, &
                      0.043108_wp, 0.043133_wp]
         case ("resp")
            charges = [ &
                      -0.139679_wp, -0.141239_wp, -0.137141_wp, -0.14177_wp, -0.139853_wp, &
                      -0.14191_wp, 0.140181_wp, 0.140495_wp, 0.138883_wp, 0.140616_wp, &
                      0.140452_wp, 0.140965_wp]
         case ("mbis")
            charges = [ &
                      -0.153068_wp, -0.15303_wp, -0.152968_wp, -0.153028_wp, -0.153022_wp, &
                      -0.152996_wp, 0.153024_wp, 0.153004_wp, 0.153033_wp, 0.153009_wp, &
                      0.153_wp, 0.153042_wp]
         case ("chelpg")
            charges = [ &
                      -0.139694_wp, -0.141476_wp, -0.137135_wp, -0.142008_wp, -0.13987_wp, &
                      -0.14213_wp, 0.140282_wp, 0.140637_wp, 0.138977_wp, 0.140759_wp, &
                      0.140555_wp, 0.141103_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (67)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.017877_wp, -0.025116_wp, -0.0131_wp, 0.01524_wp, -0.0131_wp, &
                      -0.025117_wp, 0.044709_wp, 0.044874_wp, 0.051137_wp, 0.051137_wp, &
                      0.044873_wp, 0.077519_wp, -0.235181_wp]
         case ("resp")
            charges = [ &
                      -0.108445_wp, -0.111999_wp, -0.179754_wp, 0.129902_wp, -0.181256_wp, &
                      -0.111659_wp, 0.131743_wp, 0.135448_wp, 0.160627_wp, 0.161611_wp, &
                      0.135369_wp, 0.265263_wp, -0.426848_wp]
         case ("mbis")
            charges = [ &
                      -0.10556_wp, -0.138422_wp, -0.122553_wp, -0.021603_wp, -0.122551_wp, &
                      -0.138427_wp, 0.141865_wp, 0.146821_wp, 0.157119_wp, 0.157117_wp, &
                      0.146821_wp, 0.336655_wp, -0.437284_wp]
         case ("chelpg")
            charges = [ &
                      -0.119454_wp, -0.100345_wp, -0.19572_wp, 0.148697_wp, -0.197195_wp, &
                      -0.100003_wp, 0.13364_wp, 0.133736_wp, 0.163946_wp, 0.164914_wp, &
                      0.133655_wp, 0.260019_wp, -0.42589_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.01339_wp, -0.025127_wp, -0.01197_wp, 0.006586_wp, -0.01197_wp, &
                      -0.025128_wp, 0.055987_wp, 0.055123_wp, 0.059938_wp, 0.059939_wp, &
                      0.055123_wp, 0.093665_wp, -0.298778_wp]
         case ("resp")
            charges = [ &
                      -0.106312_wp, -0.123987_wp, -0.174494_wp, 0.072592_wp, -0.176211_wp, &
                      -0.123513_wp, 0.151075_wp, 0.154959_wp, 0.176214_wp, 0.177284_wp, &
                      0.154845_wp, 0.357544_wp, -0.539995_wp]
         case ("mbis")
            charges = [ &
                      -0.102432_wp, -0.153253_wp, -0.119015_wp, -0.066012_wp, -0.119012_wp, &
                      -0.153258_wp, 0.160958_wp, 0.166985_wp, 0.172593_wp, 0.172593_wp, &
                      0.166985_wp, 0.421199_wp, -0.548333_wp]
         case ("chelpg")
            charges = [ &
                      -0.115388_wp, -0.114353_wp, -0.18781_wp, 0.088105_wp, -0.18951_wp, &
                      -0.113872_wp, 0.152638_wp, 0.153571_wp, 0.178992_wp, 0.180051_wp, &
                      0.153454_wp, 0.35345_wp, -0.539328_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.013011_wp, -0.025027_wp, -0.011925_wp, 0.0061_wp, -0.011925_wp, &
                      -0.025027_wp, 0.056591_wp, 0.055647_wp, 0.06031_wp, 0.060311_wp, &
                      0.055647_wp, 0.09443_wp, -0.302125_wp]
         case ("resp")
            charges = [ &
                      -0.105926_wp, -0.124338_wp, -0.174205_wp, 0.069411_wp, -0.175935_wp, &
                      -0.123856_wp, 0.152011_wp, 0.155873_wp, 0.176865_wp, 0.17794_wp, &
                      0.155757_wp, 0.362364_wp, -0.54596_wp]
         case ("mbis")
            charges = [ &
                      -0.101999_wp, -0.153818_wp, -0.118793_wp, -0.068434_wp, -0.11879_wp, &
                      -0.153823_wp, 0.161893_wp, 0.167953_wp, 0.173241_wp, 0.173241_wp, &
                      0.167953_wp, 0.425496_wp, -0.554123_wp]
         case ("chelpg")
            charges = [ &
                      -0.114858_wp, -0.114853_wp, -0.187316_wp, 0.084658_wp, -0.18903_wp, &
                      -0.114364_wp, 0.153551_wp, 0.154508_wp, 0.179603_wp, 0.180668_wp, &
                      0.154389_wp, 0.358365_wp, -0.545323_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (68)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.035208_wp, -0.03401_wp, -0.032583_wp, -0.005638_wp, -0.04069_wp, &
                      -0.036728_wp, 0.036767_wp, 0.037252_wp, 0.038904_wp, 0.042072_wp, &
                      0.032089_wp, 0.036341_wp, -0.247921_wp, 0.037035_wp, 0.024132_wp, &
                      0.148183_wp]
         case ("resp")
            charges = [ &
                      -0.117137_wp, -0.104794_wp, -0.210858_wp, 0.051487_wp, -0.193438_wp, &
                      -0.125367_wp, 0.120995_wp, 0.121241_wp, 0.162645_wp, 0.141067_wp, &
                      0.136791_wp, 0.126807_wp, -0.594204_wp, 0.083839_wp, 0.013402_wp, &
                      0.387525_wp]
         case ("mbis")
            charges = [ &
                      -0.137598_wp, -0.132319_wp, -0.151667_wp, -0.002175_wp, -0.174418_wp, &
                      -0.128058_wp, 0.136425_wp, 0.136476_wp, 0.14594_wp, 0.051488_wp, &
                      0.132507_wp, 0.135196_wp, -0.60906_wp, 0.106884_wp, 0.07049_wp, &
                      0.419886_wp]
         case ("chelpg")
            charges = [ &
                      -0.126213_wp, -0.095019_wp, -0.222451_wp, 0.062492_wp, -0.205684_wp, &
                      -0.115364_wp, 0.122522_wp, 0.119665_wp, 0.164838_wp, 0.143407_wp, &
                      0.139261_wp, 0.125202_wp, -0.595104_wp, 0.082705_wp, 0.012202_wp, &
                      0.387542_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.041799_wp, -0.041283_wp, -0.041538_wp, -0.008293_wp, -0.042056_wp, &
                      -0.042197_wp, 0.043329_wp, 0.043452_wp, 0.041227_wp, 0.041923_wp, &
                      0.041935_wp, 0.043365_wp, -0.272191_wp, 0.039079_wp, 0.033022_wp, &
                      0.162021_wp]
         case ("resp")
            charges = [ &
                      -0.12997_wp, -0.128102_wp, -0.229151_wp, 0.040418_wp, -0.200922_wp, &
                      -0.148424_wp, 0.137843_wp, 0.139934_wp, 0.174398_wp, 0.180168_wp, &
                      0.158345_wp, 0.145225_wp, -0.679963_wp, 0.079805_wp, 0.023134_wp, &
                      0.437263_wp]
         case ("mbis")
            charges = [ &
                      -0.154875_wp, -0.147582_wp, -0.176755_wp, 0.001424_wp, -0.184799_wp, &
                      -0.145698_wp, 0.153125_wp, 0.152554_wp, 0.156895_wp, 0.05525_wp, &
                      0.152288_wp, 0.152313_wp, -0.67183_wp, 0.11104_wp, 0.086445_wp, &
                      0.460203_wp]
         case ("chelpg")
            charges = [ &
                      -0.138071_wp, -0.119407_wp, -0.239538_wp, 0.050136_wp, -0.212125_wp, &
                      -0.139381_wp, 0.139204_wp, 0.138559_wp, 0.176361_wp, 0.183159_wp, &
                      0.160637_wp, 0.143791_wp, -0.681_wp, 0.078562_wp, 0.021809_wp, &
                      0.437304_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.042446_wp, -0.042097_wp, -0.04258_wp, -0.008563_wp, -0.042008_wp, &
                      -0.042631_wp, 0.043978_wp, 0.044009_wp, 0.041288_wp, 0.041913_wp, &
                      0.043075_wp, 0.044119_wp, -0.274656_wp, 0.039154_wp, 0.034064_wp, &
                      0.163377_wp]
         case ("resp")
            charges = [ &
                      -0.131239_wp, -0.130555_wp, -0.231054_wp, 0.038813_wp, -0.20113_wp, &
                      -0.15059_wp, 0.139514_wp, 0.14175_wp, 0.175296_wp, 0.184693_wp, &
                      0.160624_wp, 0.147082_wp, -0.688883_wp, 0.079038_wp, 0.024314_wp, &
                      0.442328_wp]
         case ("mbis")
            charges = [ &
                      -0.156578_wp, -0.149187_wp, -0.179488_wp, 0.001725_wp, -0.185478_wp, &
                      -0.14733_wp, 0.15477_wp, 0.154088_wp, 0.157767_wp, 0.055726_wp, &
                      0.154426_wp, 0.154054_wp, -0.678156_wp, 0.111205_wp, 0.088298_wp, &
                      0.464155_wp]
         case ("chelpg")
            charges = [ &
                      -0.139227_wp, -0.121987_wp, -0.241296_wp, 0.048363_wp, -0.212197_wp, &
                      -0.141665_wp, 0.140857_wp, 0.140398_wp, 0.177231_wp, 0.187775_wp, &
                      0.162894_wp, 0.145668_wp, -0.689938_wp, 0.077779_wp, 0.022971_wp, &
                      0.442372_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (69)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.026885_wp, -0.033952_wp, 0.025112_wp, -0.033952_wp, -0.026885_wp, &
                      -0.033034_wp, -0.08257_wp, 0.04161_wp, 0.044515_wp, 0.044516_wp, &
                      0.04161_wp, 0.039922_wp]
         case ("resp")
            charges = [ &
                      -0.165737_wp, -0.035566_wp, -0.01802_wp, -0.032229_wp, -0.170415_wp, &
                      -0.103732_wp, -0.107069_wp, 0.137652_wp, 0.113691_wp, 0.112903_wp, &
                      0.139563_wp, 0.128958_wp]
         case ("mbis")
            charges = [ &
                      -0.122677_wp, -0.157803_wp, 0.08836_wp, -0.157802_wp, -0.122681_wp, &
                      -0.141929_wp, -0.136383_wp, 0.140906_wp, 0.163866_wp, 0.163865_wp, &
                      0.140906_wp, 0.14138_wp]
         case ("chelpg")
            charges = [ &
                      -0.169053_wp, -0.033535_wp, -0.019023_wp, -0.030117_wp, -0.173779_wp, &
                      -0.100677_wp, -0.107189_wp, 0.138482_wp, 0.113406_wp, 0.112588_wp, &
                      0.140402_wp, 0.128497_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.028862_wp, -0.037565_wp, 0.017933_wp, -0.037565_wp, -0.028862_wp, &
                      -0.03447_wp, -0.098202_wp, 0.049567_wp, 0.050446_wp, 0.050448_wp, &
                      0.049566_wp, 0.047573_wp]
         case ("resp")
            charges = [ &
                      -0.175549_wp, -0.050097_wp, -0.02122_wp, -0.046846_wp, -0.180108_wp, &
                      -0.116229_wp, -0.127324_wp, 0.155229_wp, 0.129866_wp, 0.129108_wp, &
                      0.157086_wp, 0.146085_wp]
         case ("mbis")
            charges = [ &
                      -0.132197_wp, -0.169945_wp, 0.081884_wp, -0.169944_wp, -0.1322_wp, &
                      -0.151671_wp, -0.153316_wp, 0.157355_wp, 0.177713_wp, 0.177713_wp, &
                      0.157355_wp, 0.157261_wp]
         case ("chelpg")
            charges = [ &
                      -0.177313_wp, -0.049711_wp, -0.020737_wp, -0.046394_wp, -0.181912_wp, &
                      -0.114658_wp, -0.127599_wp, 0.155778_wp, 0.1299_wp, 0.129118_wp, &
                      0.157644_wp, 0.145882_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.029129_wp, -0.038653_wp, 0.015798_wp, -0.038653_wp, -0.02913_wp, &
                      -0.034463_wp, -0.102853_wp, 0.051791_wp, 0.05186_wp, 0.051862_wp, &
                      0.051791_wp, 0.049785_wp]
         case ("resp")
            charges = [ &
                      -0.17741_wp, -0.054815_wp, -0.021654_wp, -0.051591_wp, -0.181934_wp, &
                      -0.119059_wp, -0.133473_wp, 0.159921_wp, 0.13413_wp, 0.133382_wp, &
                      0.161763_wp, 0.15074_wp]
         case ("mbis")
            charges = [ &
                      -0.134227_wp, -0.17338_wp, 0.079785_wp, -0.173379_wp, -0.13423_wp, &
                      -0.153642_wp, -0.158298_wp, 0.161753_wp, 0.181155_wp, 0.181155_wp, &
                      0.161753_wp, 0.161563_wp]
         case ("chelpg")
            charges = [ &
                      -0.178701_wp, -0.054931_wp, -0.020712_wp, -0.051645_wp, -0.183262_wp, &
                      -0.117936_wp, -0.1338_wp, 0.160386_wp, 0.134262_wp, 0.133491_wp, &
                      0.162235_wp, 0.150613_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (70)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.064952_wp, 0.034597_wp, 0.034601_wp, 0.034507_wp, 0.01774_wp, &
                      0.042025_wp, -0.140534_wp, 0.042026_wp]
         case ("resp")
            charges = [ &
                      -0.168477_wp, 0.075607_wp, 0.076516_wp, 0.063835_wp, -0.034006_wp, &
                      0.09652_wp, -0.206948_wp, 0.096952_wp]
         case ("mbis")
            charges = [ &
                      -0.420068_wp, 0.157192_wp, 0.157193_wp, 0.14224_wp, -0.139055_wp, &
                      0.153466_wp, -0.204418_wp, 0.153461_wp]
         case ("chelpg")
            charges = [ &
                      -0.176488_wp, 0.077819_wp, 0.078715_wp, 0.065896_wp, -0.032842_wp, &
                      0.096712_wp, -0.206965_wp, 0.097152_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.063872_wp, 0.035383_wp, 0.035388_wp, 0.043335_wp, 0.026437_wp, &
                      0.050992_wp, -0.178647_wp, 0.050994_wp]
         case ("resp")
            charges = [ &
                      -0.189693_wp, 0.078642_wp, 0.079524_wp, 0.080027_wp, 0.020504_wp, &
                      0.101137_wp, -0.271794_wp, 0.101653_wp]
         case ("mbis")
            charges = [ &
                      -0.434311_wp, 0.159675_wp, 0.159678_wp, 0.158586_wp, -0.119859_wp, &
                      0.165478_wp, -0.254711_wp, 0.165475_wp]
         case ("chelpg")
            charges = [ &
                      -0.199751_wp, 0.081252_wp, 0.082117_wp, 0.082512_wp, 0.024099_wp, &
                      0.100683_wp, -0.272122_wp, 0.10121_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.063616_wp, 0.035401_wp, 0.035406_wp, 0.044898_wp, 0.027957_wp, &
                      0.052448_wp, -0.184932_wp, 0.05245_wp]
         case ("resp")
            charges = [ &
                      -0.193897_wp, 0.079111_wp, 0.079985_wp, 0.083047_wp, 0.030671_wp, &
                      0.101684_wp, -0.282816_wp, 0.102215_wp]
         case ("mbis")
            charges = [ &
                      -0.436221_wp, 0.159805_wp, 0.159808_wp, 0.161395_wp, -0.116453_wp, &
                      0.167356_wp, -0.263032_wp, 0.167353_wp]
         case ("chelpg")
            charges = [ &
                      -0.204308_wp, 0.08179_wp, 0.082646_wp, 0.085603_wp, 0.034696_wp, &
                      0.101114_wp, -0.2832_wp, 0.101658_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (71)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.075717_wp, -0.045252_wp, -0.045285_wp, -0.045141_wp, 0.059962_wp]
         case ("resp")
            charges = [ &
                      -0.365361_wp, 0.027781_wp, 0.028633_wp, 0.028385_wp, 0.280562_wp]
         case ("mbis")
            charges = [ &
                      -0.1375_wp, -0.031573_wp, -0.031634_wp, -0.03141_wp, 0.232117_wp]
         case ("chelpg")
            charges = [ &
                      -0.368243_wp, 0.028393_wp, 0.029239_wp, 0.028991_wp, 0.281621_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.079978_wp, -0.05133_wp, -0.05136_wp, -0.051181_wp, 0.073893_wp]
         case ("resp")
            charges = [ &
                      -0.392424_wp, 0.021551_wp, 0.022419_wp, 0.022294_wp, 0.32616_wp]
         case ("mbis")
            charges = [ &
                      -0.128618_wp, -0.041661_wp, -0.041719_wp, -0.041449_wp, 0.253448_wp]
         case ("chelpg")
            charges = [ &
                      -0.395254_wp, 0.022151_wp, 0.023013_wp, 0.02289_wp, 0.327199_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.081789_wp, -0.053701_wp, -0.05373_wp, -0.053533_wp, 0.079176_wp]
         case ("resp")
            charges = [ &
                      -0.403403_wp, 0.018809_wp, 0.019681_wp, 0.019618_wp, 0.345295_wp]
         case ("mbis")
            charges = [ &
                      -0.124549_wp, -0.045644_wp, -0.0457_wp, -0.045406_wp, 0.261301_wp]
         case ("chelpg")
            charges = [ &
                      -0.40621_wp, 0.019402_wp, 0.020269_wp, 0.020206_wp, 0.346333_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (72)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069042_wp, -0.031727_wp, -0.037762_wp, 0.157709_wp, -0.254809_wp, &
                      0.02253_wp, 0.028006_wp, 0.022545_wp, 0.027691_wp, 0.027703_wp, &
                      0.038472_wp, 0.038398_wp, 0.03028_wp]
         case ("resp")
            charges = [ &
                      -0.423659_wp, 0.287095_wp, -0.140871_wp, 0.433515_wp, -0.474586_wp, &
                      0.09792_wp, 0.097964_wp, 0.097661_wp, -0.022103_wp, -0.021846_wp, &
                      0.040354_wp, 0.040364_wp, -0.011808_wp]
         case ("mbis")
            charges = [ &
                      -0.428504_wp, -0.126271_wp, -0.387198_wp, 0.462723_wp, -0.461824_wp, &
                      0.128868_wp, 0.13701_wp, 0.128867_wp, 0.115794_wp, 0.115929_wp, &
                      0.148489_wp, 0.148357_wp, 0.017753_wp]
         case ("chelpg")
            charges = [ &
                      -0.44286_wp, 0.305401_wp, -0.159521_wp, 0.441424_wp, -0.476743_wp, &
                      0.102495_wp, 0.101894_wp, 0.102243_wp, -0.024948_wp, -0.024686_wp, &
                      0.044236_wp, 0.044245_wp, -0.013182_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069484_wp, -0.033784_wp, -0.03371_wp, 0.169715_wp, -0.303846_wp, &
                      0.025412_wp, 0.027409_wp, 0.025402_wp, 0.023821_wp, 0.023775_wp, &
                      0.050053_wp, 0.050003_wp, 0.045227_wp]
         case ("resp")
            charges = [ &
                      -0.429539_wp, 0.309376_wp, -0.139269_wp, 0.476837_wp, -0.563823_wp, &
                      0.101855_wp, 0.095504_wp, 0.101565_wp, -0.035357_wp, -0.035235_wp, &
                      0.055169_wp, 0.055239_wp, 0.007678_wp]
         case ("mbis")
            charges = [ &
                      -0.43766_wp, -0.11968_wp, -0.409895_wp, 0.503383_wp, -0.539676_wp, &
                      0.135178_wp, 0.137099_wp, 0.135132_wp, 0.10892_wp, 0.108957_wp, &
                      0.169782_wp, 0.169694_wp, 0.038759_wp]
         case ("chelpg")
            charges = [ &
                      -0.449039_wp, 0.328271_wp, -0.158355_wp, 0.484798_wp, -0.56598_wp, &
                      0.106482_wp, 0.09948_wp, 0.1062_wp, -0.038328_wp, -0.038202_wp, &
                      0.059141_wp, 0.059209_wp, 0.006323_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069469_wp, -0.034027_wp, -0.033252_wp, 0.170885_wp, -0.308858_wp, &
                      0.025789_wp, 0.027244_wp, 0.025777_wp, 0.023263_wp, 0.02321_wp, &
                      0.05134_wp, 0.051293_wp, 0.046797_wp]
         case ("resp")
            charges = [ &
                      -0.429194_wp, 0.310988_wp, -0.138903_wp, 0.480972_wp, -0.572885_wp, &
                      0.10225_wp, 0.094845_wp, 0.101955_wp, -0.0369_wp, -0.036796_wp, &
                      0.056886_wp, 0.056963_wp, 0.009816_wp]
         case ("mbis")
            charges = [ &
                      -0.43838_wp, -0.118845_wp, -0.412173_wp, 0.50736_wp, -0.547626_wp, &
                      0.135928_wp, 0.136866_wp, 0.135876_wp, 0.107914_wp, 0.107938_wp, &
                      0.17211_wp, 0.172028_wp, 0.040998_wp]
         case ("chelpg")
            charges = [ &
                      -0.44872_wp, 0.329941_wp, -0.158029_wp, 0.488936_wp, -0.575041_wp, &
                      0.106881_wp, 0.098824_wp, 0.106595_wp, -0.039884_wp, -0.039776_wp, &
                      0.060866_wp, 0.060942_wp, 0.008464_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (73)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.23362_wp, -0.299246_wp, -0.190842_wp, 0.178519_wp, -0.034245_wp, &
                      -0.030729_wp, 0.040576_wp, 0.040847_wp, -0.06796_wp, 0.027523_wp, &
                      0.027481_wp, 0.02827_wp, 0.023085_wp, 0.023099_wp]
         case ("resp")
            charges = [ &
                      0.6542_wp, -0.557141_wp, -0.607222_wp, 0.431337_wp, -0.165554_wp, &
                      0.163073_wp, 0.08386_wp, 0.082917_wp, -0.392199_wp, 0.009922_wp, &
                      0.00895_wp, 0.103017_wp, 0.092544_wp, 0.092298_wp]
         case ("mbis")
            charges = [ &
                      0.8305_wp, -0.601086_wp, -0.653635_wp, 0.459901_wp, -0.446026_wp, &
                      -0.125803_wp, 0.167219_wp, 0.167071_wp, -0.42644_wp, 0.116372_wp, &
                      0.116071_wp, 0.137118_wp, 0.12934_wp, 0.129397_wp]
         case ("chelpg")
            charges = [ &
                      0.66252_wp, -0.559273_wp, -0.608_wp, 0.430952_wp, -0.186405_wp, &
                      0.1805_wp, 0.088404_wp, 0.087458_wp, -0.410594_wp, 0.00745_wp, &
                      0.006455_wp, 0.106723_wp, 0.097014_wp, 0.096796_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.238181_wp, -0.324808_wp, -0.191981_wp, 0.189952_wp, -0.033632_wp, &
                      -0.031544_wp, 0.045802_wp, 0.046169_wp, -0.068184_wp, 0.025745_wp, &
                      0.025806_wp, 0.028493_wp, 0.024994_wp, 0.025005_wp]
         case ("resp")
            charges = [ &
                      0.701985_wp, -0.609607_wp, -0.635637_wp, 0.460419_wp, -0.193956_wp, &
                      0.201026_wp, 0.094788_wp, 0.093904_wp, -0.392605_wp, -0.003684_wp, &
                      -0.004469_wp, 0.09946_wp, 0.094346_wp, 0.094029_wp]
         case ("mbis")
            charges = [ &
                      0.863386_wp, -0.643367_wp, -0.674539_wp, 0.48464_wp, -0.466041_wp, &
                      -0.11916_wp, 0.178665_wp, 0.178685_wp, -0.433036_wp, 0.112793_wp, &
                      0.112681_wp, 0.1382_wp, 0.133519_wp, 0.133572_wp]
         case ("chelpg")
            charges = [ &
                      0.710543_wp, -0.611774_wp, -0.636367_wp, 0.459981_wp, -0.215656_wp, &
                      0.219495_wp, 0.099483_wp, 0.098598_wp, -0.411431_wp, -0.006386_wp, &
                      -0.007195_wp, 0.103215_wp, 0.09889_wp, 0.098603_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.240676_wp, -0.340407_wp, -0.191583_wp, 0.196201_wp, -0.033115_wp, &
                      -0.032134_wp, 0.049225_wp, 0.049644_wp, -0.068098_wp, 0.024117_wp, &
                      0.024249_wp, 0.028345_wp, 0.026434_wp, 0.026444_wp]
         case ("resp")
            charges = [ &
                      0.728537_wp, -0.641142_wp, -0.649434_wp, 0.476216_wp, -0.209386_wp, &
                      0.221856_wp, 0.101645_wp, 0.100789_wp, -0.389621_wp, -0.012608_wp, &
                      -0.013262_wp, 0.096078_wp, 0.095347_wp, 0.094986_wp]
         case ("mbis")
            charges = [ &
                      0.881721_wp, -0.668819_wp, -0.683927_wp, 0.497516_wp, -0.4777_wp, &
                      -0.114631_wp, 0.185898_wp, 0.186008_wp, -0.436299_wp, 0.109608_wp, &
                      0.109625_wp, 0.138154_wp, 0.136397_wp, 0.136448_wp]
         case ("chelpg")
            charges = [ &
                      0.737211_wp, -0.643322_wp, -0.650149_wp, 0.475757_wp, -0.231515_wp, &
                      0.240875_wp, 0.106416_wp, 0.105559_wp, -0.408638_wp, -0.015439_wp, &
                      -0.016118_wp, 0.099847_wp, 0.099922_wp, 0.099592_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (74)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068829_wp, -0.039038_wp, 0.179186_wp, -0.270294_wp, -0.07293_wp, &
                      0.029061_wp, 0.024837_wp, 0.029028_wp, 0.034859_wp, 0.034748_wp, &
                      0.038982_wp, 0.039096_wp, 0.041292_wp]
         case ("resp")
            charges = [ &
                      -0.305665_wp, 0.005992_wp, 0.660217_wp, -0.525239_wp, -0.616866_wp, &
                      0.087697_wp, 0.081416_wp, 0.086915_wp, 0.022304_wp, 0.021805_wp, &
                      0.148258_wp, 0.149174_wp, 0.183993_wp]
         case ("mbis")
            charges = [ &
                      -0.389183_wp, -0.355947_wp, 0.649229_wp, -0.506235_wp, -0.631065_wp, &
                      0.14253_wp, 0.126475_wp, 0.142393_wp, 0.141107_wp, 0.14106_wp, &
                      0.173924_wp, 0.174064_wp, 0.191647_wp]
         case ("chelpg")
            charges = [ &
                      -0.317508_wp, 0.005957_wp, 0.66837_wp, -0.526597_wp, -0.636126_wp, &
                      0.090613_wp, 0.084858_wp, 0.089841_wp, 0.022837_wp, 0.022332_wp, &
                      0.152717_wp, 0.153642_wp, 0.189063_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072462_wp, -0.035206_wp, 0.188383_wp, -0.321142_wp, -0.071038_wp, &
                      0.024194_wp, 0.03031_wp, 0.024227_wp, 0.046756_wp, 0.046724_wp, &
                      0.050374_wp, 0.050522_wp, 0.038357_wp]
         case ("resp")
            charges = [ &
                      -0.324194_wp, 0.028449_wp, 0.703902_wp, -0.616476_wp, -0.642124_wp, &
                      0.080927_wp, 0.095606_wp, 0.080192_wp, 0.034404_wp, 0.034082_wp, &
                      0.170707_wp, 0.171674_wp, 0.182851_wp]
         case ("mbis")
            charges = [ &
                      -0.395973_wp, -0.379012_wp, 0.703441_wp, -0.587794_wp, -0.661301_wp, &
                      0.135993_wp, 0.138418_wp, 0.135976_wp, 0.163441_wp, 0.163545_wp, &
                      0.196641_wp, 0.196837_wp, 0.189785_wp]
         case ("chelpg")
            charges = [ &
                      -0.336941_wp, 0.029655_wp, 0.711612_wp, -0.617762_wp, -0.661288_wp, &
                      0.084016_wp, 0.09922_wp, 0.083291_wp, 0.03466_wp, 0.034332_wp, &
                      0.175146_wp, 0.176123_wp, 0.187936_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072715_wp, -0.034877_wp, 0.188975_wp, -0.324943_wp, -0.070862_wp, &
                      0.02371_wp, 0.0308_wp, 0.02375_wp, 0.047725_wp, 0.047703_wp, &
                      0.051289_wp, 0.051439_wp, 0.038008_wp]
         case ("resp")
            charges = [ &
                      -0.325742_wp, 0.030783_wp, 0.70618_wp, -0.623131_wp, -0.643202_wp, &
                      0.080227_wp, 0.096881_wp, 0.079498_wp, 0.035361_wp, 0.035058_wp, &
                      0.172373_wp, 0.173344_wp, 0.18237_wp]
         case ("mbis")
            charges = [ &
                      -0.39621_wp, -0.380672_wp, 0.707173_wp, -0.593851_wp, -0.663269_wp, &
                      0.135245_wp, 0.139422_wp, 0.13524_wp, 0.165219_wp, 0.165339_wp, &
                      0.198398_wp, 0.198598_wp, 0.189366_wp]
         case ("chelpg")
            charges = [ &
                      -0.338575_wp, 0.032116_wp, 0.713842_wp, -0.62441_wp, -0.662348_wp, &
                      0.083332_wp, 0.10051_wp, 0.082614_wp, 0.035589_wp, 0.035279_wp, &
                      0.176808_wp, 0.177789_wp, 0.187454_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (75)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.249993_wp, 0.081274_wp, -0.002693_wp, 0.052222_wp, 0.052225_wp, &
                      -0.021189_wp, -0.065035_wp, 0.033239_wp, 0.033241_wp, 0.032403_wp, &
                      0.027155_wp, 0.027152_wp]
         case ("resp")
            charges = [ &
                      -0.452193_wp, 0.281949_wp, -0.15466_wp, 0.101738_wp, 0.101244_wp, &
                      0.187563_wp, -0.380389_wp, 0.010881_wp, 0.011032_wp, 0.102201_wp, &
                      0.095534_wp, 0.095099_wp]
         case ("mbis")
            charges = [ &
                      -0.463171_wp, 0.363904_wp, -0.314552_wp, 0.166161_wp, 0.166146_wp, &
                      -0.152139_wp, -0.425293_wp, 0.123266_wp, 0.123292_wp, 0.143192_wp, &
                      0.134594_wp, 0.1346_wp]
         case ("chelpg")
            charges = [ &
                      -0.454667_wp, 0.29052_wp, -0.177842_wp, 0.106909_wp, 0.106446_wp, &
                      0.204498_wp, -0.395541_wp, 0.008144_wp, 0.00832_wp, 0.105098_wp, &
                      0.09927_wp, 0.098847_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.320025_wp, 0.102984_wp, 0.003441_wp, 0.066391_wp, 0.066395_wp, &
                      -0.018268_wp, -0.064166_wp, 0.034939_wp, 0.03493_wp, 0.032845_wp, &
                      0.030269_wp, 0.030267_wp]
         case ("resp")
            charges = [ &
                      -0.570166_wp, 0.36247_wp, -0.199432_wp, 0.129724_wp, 0.129329_wp, &
                      0.237443_wp, -0.39027_wp, -1.6e-05_wp, 0.000209_wp, 0.098917_wp, &
                      0.101107_wp, 0.100684_wp]
         case ("mbis")
            charges = [ &
                      -0.586353_wp, 0.456642_wp, -0.355943_wp, 0.193545_wp, 0.193534_wp, &
                      -0.141554_wp, -0.433665_wp, 0.123964_wp, 0.123968_wp, 0.144433_wp, &
                      0.14071_wp, 0.140719_wp]
         case ("chelpg")
            charges = [ &
                      -0.572782_wp, 0.371583_wp, -0.224082_wp, 0.135208_wp, 0.134847_wp, &
                      0.255426_wp, -0.405759_wp, -0.002967_wp, -0.002716_wp, 0.101826_wp, &
                      0.104913_wp, 0.104503_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.323868_wp, 0.104171_wp, 0.0038_wp, 0.06719_wp, 0.067195_wp, &
                      -0.018113_wp, -0.064095_wp, 0.034987_wp, 0.034977_wp, 0.032824_wp, &
                      0.030467_wp, 0.030465_wp]
         case ("resp")
            charges = [ &
                      -0.576717_wp, 0.366988_wp, -0.201881_wp, 0.131312_wp, 0.130924_wp, &
                      0.239967_wp, -0.390406_wp, -0.000665_wp, -0.000437_wp, 0.098555_wp, &
                      0.10139_wp, 0.100968_wp]
         case ("mbis")
            charges = [ &
                      -0.593103_wp, 0.461707_wp, -0.358172_wp, 0.195068_wp, 0.195056_wp, &
                      -0.140931_wp, -0.433998_wp, 0.123917_wp, 0.12392_wp, 0.144389_wp, &
                      0.141069_wp, 0.141077_wp]
         case ("chelpg")
            charges = [ &
                      -0.579348_wp, 0.376126_wp, -0.226595_wp, 0.13681_wp, 0.136455_wp, &
                      0.258016_wp, -0.405903_wp, -0.003631_wp, -0.003376_wp, 0.101463_wp, &
                      0.105196_wp, 0.104787_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (76)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.303487_wp, 0.229963_wp, -0.148224_wp, 0.052155_wp, -0.037858_wp, &
                      -0.030471_wp, -0.067664_wp, -0.066045_wp, 0.032959_wp, 0.032641_wp, &
                      0.025933_wp, 0.026072_wp, 0.023107_wp, 0.023035_wp, 0.023879_wp, &
                      0.023918_wp, 0.026091_wp, 0.04782_wp, 0.04568_wp, 0.040494_wp]
         case ("resp")
            charges = [ &
                      -0.601184_wp, 0.912293_wp, -0.469783_wp, 0.124238_wp, -0.076783_wp, &
                      0.142739_wp, -0.294935_wp, -0.692004_wp, 0.051657_wp, 0.050343_wp, &
                      0.050967_wp, 0.050819_wp, -0.014741_wp, -0.014349_wp, 0.067792_wp, &
                      0.068341_wp, 0.078886_wp, 0.190017_wp, 0.189695_wp, 0.185994_wp]
         case ("mbis")
            charges = [ &
                      -0.601597_wp, 0.843137_wp, -0.41268_wp, 0.063187_wp, -0.290665_wp, &
                      -0.14979_wp, -0.420971_wp, -0.666489_wp, 0.096675_wp, 0.096065_wp, &
                      0.130199_wp, 0.130393_wp, 0.103744_wp, 0.103643_wp, 0.13052_wp, &
                      0.130635_wp, 0.133202_wp, 0.192798_wp, 0.198024_wp, 0.189969_wp]
         case ("chelpg")
            charges = [ &
                      -0.602975_wp, 0.919976_wp, -0.472708_wp, 0.134249_wp, -0.090617_wp, &
                      0.156491_wp, -0.310396_wp, -0.709219_wp, 0.049247_wp, 0.047917_wp, &
                      0.053872_wp, 0.053599_wp, -0.016975_wp, -0.016606_wp, 0.071426_wp, &
                      0.072008_wp, 0.082209_wp, 0.194062_wp, 0.193921_wp, 0.190518_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.336445_wp, 0.234692_wp, -0.148641_wp, 0.051493_wp, -0.037555_wp, &
                      -0.029349_wp, -0.06933_wp, -0.065356_wp, 0.034841_wp, 0.034763_wp, &
                      0.027314_wp, 0.027416_wp, 0.025237_wp, 0.02521_wp, 0.024462_wp, &
                      0.024484_wp, 0.026354_wp, 0.055471_wp, 0.04851_wp, 0.046428_wp]
         case ("resp")
            charges = [ &
                      -0.668314_wp, 0.96474_wp, -0.500612_wp, 0.146301_wp, -0.070584_wp, &
                      0.164503_wp, -0.311465_wp, -0.716733_wp, 0.051221_wp, 0.050464_wp, &
                      0.049613_wp, 0.049209_wp, -0.018052_wp, -0.017644_wp, 0.070709_wp, &
                      0.071241_wp, 0.080868_wp, 0.205742_wp, 0.198057_wp, 0.200736_wp]
         case ("mbis")
            charges = [ &
                      -0.654859_wp, 0.87926_wp, -0.420541_wp, 0.058638_wp, -0.292207_wp, &
                      -0.147046_wp, -0.432837_wp, -0.692063_wp, 0.101141_wp, 0.100983_wp, &
                      0.132718_wp, 0.132847_wp, 0.106843_wp, 0.106823_wp, 0.133588_wp, &
                      0.133673_wp, 0.1357_wp, 0.208672_wp, 0.205712_wp, 0.202956_wp]
         case ("chelpg")
            charges = [ &
                      -0.670101_wp, 0.972447_wp, -0.503608_wp, 0.156643_wp, -0.084463_wp, &
                      0.178589_wp, -0.327168_wp, -0.734024_wp, 0.048708_wp, 0.047937_wp, &
                      0.052483_wp, 0.051953_wp, -0.020379_wp, -0.019993_wp, 0.074392_wp, &
                      0.074957_wp, 0.08424_wp, 0.209805_wp, 0.202299_wp, 0.205284_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.346436_wp, 0.235914_wp, -0.148364_wp, 0.051306_wp, -0.037356_wp, &
                      -0.028972_wp, -0.069748_wp, -0.065034_wp, 0.035182_wp, 0.035204_wp, &
                      0.027789_wp, 0.027859_wp, 0.025847_wp, 0.025846_wp, 0.02462_wp, &
                      0.02463_wp, 0.026359_wp, 0.057804_wp, 0.049153_wp, 0.048396_wp]
         case ("resp")
            charges = [ &
                      -0.688349_wp, 0.978568_wp, -0.508269_wp, 0.152563_wp, -0.068308_wp, &
                      0.170936_wp, -0.315852_wp, -0.722591_wp, 0.050737_wp, 0.050255_wp, &
                      0.049244_wp, 0.048688_wp, -0.019065_wp, -0.01865_wp, 0.071443_wp, &
                      0.071961_wp, 0.081225_wp, 0.210337_wp, 0.199901_wp, 0.205224_wp]
         case ("mbis")
            charges = [ &
                      -0.67086_wp, 0.889258_wp, -0.422062_wp, 0.057531_wp, -0.29239_wp, &
                      -0.146241_wp, -0.435891_wp, -0.698981_wp, 0.102045_wp, 0.102073_wp, &
                      0.133504_wp, 0.133575_wp, 0.107708_wp, 0.107736_wp, 0.134392_wp, &
                      0.134453_wp, 0.136227_wp, 0.213361_wp, 0.20753_wp, 0.207032_wp]
         case ("chelpg")
            charges = [ &
                      -0.690136_wp, 0.986284_wp, -0.511286_wp, 0.162966_wp, -0.082153_wp, &
                      0.185096_wp, -0.331611_wp, -0.739912_wp, 0.048202_wp, 0.047706_wp, &
                      0.052095_wp, 0.051412_wp, -0.021414_wp, -0.021022_wp, 0.075136_wp, &
                      0.075686_wp, 0.084609_wp, 0.214409_wp, 0.20415_wp, 0.209783_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (77)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.009737_wp, -0.044365_wp, -0.233401_wp, 0.016426_wp, 0.016417_wp, &
                      -0.033442_wp, 0.024103_wp, 0.024059_wp, 0.090504_wp, 0.090513_wp, &
                      -0.069643_wp, 0.020566_wp, 0.020558_wp, 0.024178_wp, 0.021894_wp, &
                      0.021899_wp]
         case ("resp")
            charges = [ &
                      0.397472_wp, -0.201845_wp, -0.91509_wp, 0.028672_wp, 0.027993_wp, &
                      0.244667_wp, -0.017996_wp, -0.018842_wp, 0.322087_wp, 0.322492_wp, &
                      -0.23129_wp, -0.053978_wp, -0.054292_wp, 0.049936_wp, 0.050291_wp, &
                      0.049723_wp]
         case ("mbis")
            charges = [ &
                      0.061963_wp, -0.316726_wp, -0.830347_wp, 0.114643_wp, 0.114617_wp, &
                      -0.136925_wp, 0.083915_wp, 0.083872_wp, 0.331044_wp, 0.331046_wp, &
                      -0.42133_wp, 0.098763_wp, 0.098738_wp, 0.130235_wp, 0.128237_wp, &
                      0.128259_wp]
         case ("chelpg")
            charges = [ &
                      0.411364_wp, -0.221131_wp, -0.919749_wp, 0.032087_wp, 0.031431_wp, &
                      0.263651_wp, -0.020783_wp, -0.021632_wp, 0.32316_wp, 0.323558_wp, &
                      -0.246686_wp, -0.057445_wp, -0.057772_wp, 0.052777_wp, 0.053868_wp, &
                      0.053302_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.008056_wp, -0.042813_wp, -0.253475_wp, 0.01939_wp, 0.019395_wp, &
                      -0.033656_wp, 0.023773_wp, 0.023712_wp, 0.097731_wp, 0.09774_wp, &
                      -0.071624_wp, 0.021353_wp, 0.021348_wp, 0.024329_wp, 0.022372_wp, &
                      0.022372_wp]
         case ("resp")
            charges = [ &
                      0.447936_wp, -0.177784_wp, -1.043826_wp, 0.028413_wp, 0.027727_wp, &
                      0.229179_wp, -0.027801_wp, -0.028816_wp, 0.365862_wp, 0.366149_wp, &
                      -0.235497_wp, -0.051578_wp, -0.051853_wp, 0.052174_wp, 0.050138_wp, &
                      0.049577_wp]
         case ("mbis")
            charges = [ &
                      0.078806_wp, -0.314359_wp, -0.914029_wp, 0.118924_wp, 0.118927_wp, &
                      -0.14058_wp, 0.082665_wp, 0.082591_wp, 0.361447_wp, 0.36145_wp, &
                      -0.432291_wp, 0.100732_wp, 0.100712_wp, 0.132719_wp, 0.131139_wp, &
                      0.131152_wp]
         case ("chelpg")
            charges = [ &
                      0.462064_wp, -0.196925_wp, -1.048523_wp, 0.031776_wp, 0.031112_wp, &
                      0.248006_wp, -0.030668_wp, -0.031688_wp, 0.366927_wp, 0.367206_wp, &
                      -0.250875_wp, -0.055011_wp, -0.055299_wp, 0.055019_wp, 0.05372_wp, &
                      0.05316_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.007478_wp, -0.042094_wp, -0.260161_wp, 0.020606_wp, 0.020615_wp, &
                      -0.033713_wp, 0.023429_wp, 0.02336_wp, 0.100137_wp, 0.100146_wp, &
                      -0.072186_wp, 0.021516_wp, 0.021511_wp, 0.024255_wp, 0.022553_wp, &
                      0.022552_wp]
         case ("resp")
            charges = [ &
                      0.464395_wp, -0.167856_wp, -1.087946_wp, 0.028495_wp, 0.027808_wp, &
                      0.222342_wp, -0.031378_wp, -0.032453_wp, 0.380921_wp, 0.381167_wp, &
                      -0.235885_wp, -0.050617_wp, -0.050878_wp, 0.052577_wp, 0.049933_wp, &
                      0.049374_wp]
         case ("mbis")
            charges = [ &
                      0.084824_wp, -0.313258_wp, -0.941783_wp, 0.120661_wp, 0.120674_wp, &
                      -0.141785_wp, 0.081797_wp, 0.08171_wp, 0.371443_wp, 0.371445_wp, &
                      -0.435401_wp, 0.10119_wp, 0.101171_wp, 0.133221_wp, 0.132042_wp, &
                      0.132053_wp]
         case ("chelpg")
            charges = [ &
                      0.478607_wp, -0.186921_wp, -1.09267_wp, 0.031834_wp, 0.03117_wp, &
                      0.241096_wp, -0.034276_wp, -0.035356_wp, 0.381988_wp, 0.382225_wp, &
                      -0.251218_wp, -0.054037_wp, -0.05431_wp, 0.055415_wp, 0.053505_wp, &
                      0.052948_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (78)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037165_wp, -0.041367_wp, 0.008482_wp, -0.041264_wp, -0.03715_wp, &
                      -0.041839_wp, -0.035799_wp, -0.035206_wp, -0.034111_wp, -0.069557_wp, &
                      0.034718_wp, 0.03174_wp, 0.031816_wp, 0.034717_wp, 0.034224_wp, &
                      0.024855_wp, 0.024865_wp, 0.019822_wp, 0.019881_wp, 0.019712_wp, &
                      0.019682_wp, 0.023963_wp, 0.022492_wp, 0.022496_wp]
         case ("resp")
            charges = [ &
                      -0.13684_wp, -0.254547_wp, 0.277956_wp, -0.254322_wp, -0.136865_wp, &
                      -0.129164_wp, -0.357735_wp, 0.141927_wp, 0.159701_wp, -0.373434_wp, &
                      0.13308_wp, 0.152069_wp, 0.152623_wp, 0.132777_wp, 0.124287_wp, &
                      0.084947_wp, 0.085676_wp, -0.010531_wp, -0.010068_wp, -0.017592_wp, &
                      -0.018825_wp, 0.092273_wp, 0.081506_wp, 0.081099_wp]
         case ("mbis")
            charges = [ &
                      -0.120159_wp, -0.200318_wp, 0.118466_wp, -0.200282_wp, -0.1202_wp, &
                      -0.155628_wp, -0.267382_wp, -0.205191_wp, -0.164653_wp, -0.420965_wp, &
                      0.13229_wp, 0.136481_wp, 0.136673_wp, 0.132294_wp, 0.135954_wp, &
                      0.122223_wp, 0.122313_wp, 0.113711_wp, 0.11379_wp, 0.101159_wp, &
                      0.101089_wp, 0.130485_wp, 0.128923_wp, 0.128934_wp]
         case ("chelpg")
            charges = [ &
                      -0.122189_wp, -0.278665_wp, 0.318057_wp, -0.278657_wp, -0.122118_wp, &
                      -0.143731_wp, -0.398055_wp, 0.15718_wp, 0.169891_wp, -0.3954_wp, &
                      0.131218_wp, 0.157069_wp, 0.157697_wp, 0.130901_wp, 0.126944_wp, &
                      0.092378_wp, 0.093132_wp, -0.013341_wp, -0.012907_wp, -0.019027_wp, &
                      -0.020296_wp, 0.097792_wp, 0.086266_wp, 0.085861_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.041968_wp, -0.045224_wp, 0.005647_wp, -0.045095_wp, -0.041954_wp, &
                      -0.047006_wp, -0.036208_wp, -0.035374_wp, -0.034_wp, -0.071017_wp, &
                      0.038477_wp, 0.036024_wp, 0.03611_wp, 0.038479_wp, 0.037775_wp, &
                      0.027056_wp, 0.02704_wp, 0.019905_wp, 0.019911_wp, 0.020801_wp, &
                      0.020787_wp, 0.024316_wp, 0.022763_wp, 0.022762_wp]
         case ("resp")
            charges = [ &
                      -0.151078_wp, -0.26833_wp, 0.278213_wp, -0.268296_wp, -0.150993_wp, &
                      -0.138869_wp, -0.371234_wp, 0.160985_wp, 0.171212_wp, -0.390316_wp, &
                      0.144357_wp, 0.16552_wp, 0.166183_wp, 0.14403_wp, 0.134375_wp, &
                      0.089916_wp, 0.090605_wp, -0.015552_wp, -0.015168_wp, -0.019213_wp, &
                      -0.020424_wp, 0.095879_wp, 0.084318_wp, 0.083881_wp]
         case ("mbis")
            charges = [ &
                      -0.131133_wp, -0.212634_wp, 0.118747_wp, -0.212529_wp, -0.131203_wp, &
                      -0.167881_wp, -0.273527_wp, -0.204988_wp, -0.164364_wp, -0.429627_wp, &
                      0.142395_wp, 0.147463_wp, 0.14766_wp, 0.142407_wp, 0.145955_wp, &
                      0.127195_wp, 0.127233_wp, 0.114039_wp, 0.114014_wp, 0.103064_wp, &
                      0.103023_wp, 0.132716_wp, 0.130989_wp, 0.13099_wp]
         case ("chelpg")
            charges = [ &
                      -0.13658_wp, -0.292389_wp, 0.318482_wp, -0.292578_wp, -0.136394_wp, &
                      -0.153317_wp, -0.412059_wp, 0.176737_wp, 0.181487_wp, -0.412481_wp, &
                      0.142536_wp, 0.170511_wp, 0.171251_wp, 0.142192_wp, 0.137017_wp, &
                      0.097435_wp, 0.098149_wp, -0.01849_wp, -0.018136_wp, -0.020686_wp, &
                      -0.021933_wp, 0.101453_wp, 0.089113_wp, 0.088679_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.046137_wp, -0.048472_wp, 0.003152_wp, -0.048316_wp, -0.046122_wp, &
                      -0.05158_wp, -0.036449_wp, -0.035452_wp, -0.033712_wp, -0.072159_wp, &
                      0.041593_wp, 0.039783_wp, 0.039881_wp, 0.041599_wp, 0.040637_wp, &
                      0.029095_wp, 0.029058_wp, 0.01977_wp, 0.019715_wp, 0.021919_wp, &
                      0.021924_wp, 0.024592_wp, 0.022851_wp, 0.022839_wp]
         case ("resp")
            charges = [ &
                      -0.16355_wp, -0.279126_wp, 0.276007_wp, -0.279246_wp, -0.163375_wp, &
                      -0.147193_wp, -0.380555_wp, 0.177092_wp, 0.182245_wp, -0.404414_wp, &
                      0.153886_wp, 0.177036_wp, 0.177799_wp, 0.15354_wp, 0.142763_wp, &
                      0.094045_wp, 0.094699_wp, -0.02037_wp, -0.020085_wp, -0.02058_wp, &
                      -0.021757_wp, 0.098828_wp, 0.086392_wp, 0.085919_wp]
         case ("mbis")
            charges = [ &
                      -0.140483_wp, -0.222909_wp, 0.118452_wp, -0.222741_wp, -0.140576_wp, &
                      -0.178578_wp, -0.278576_wp, -0.20428_wp, -0.164004_wp, -0.436323_wp, &
                      0.150849_wp, 0.156943_wp, 0.157148_wp, 0.150869_wp, 0.154257_wp, &
                      0.131657_wp, 0.131651_wp, 0.113858_wp, 0.113716_wp, 0.104939_wp, &
                      0.104934_wp, 0.13446_wp, 0.13238_wp, 0.132362_wp]
         case ("chelpg")
            charges = [ &
                      -0.14917_wp, -0.303134_wp, 0.316392_wp, -0.30348_wp, -0.148891_wp, &
                      -0.161552_wp, -0.421746_wp, 0.193221_wp, 0.192601_wp, -0.426732_wp, &
                      0.152095_wp, 0.182018_wp, 0.182859_wp, 0.151732_wp, 0.145396_wp, &
                      0.101626_wp, 0.102306_wp, -0.023408_wp, -0.023153_wp, -0.022086_wp, &
                      -0.0233_wp, 0.104444_wp, 0.091216_wp, 0.090745_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (79)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073725_wp, -0.035031_wp, -0.00339_wp, -0.071086_wp, 0.010961_wp, &
                      -0.039811_wp, -0.036615_wp, -0.041775_wp, -0.036318_wp, -0.0375_wp, &
                      0.021029_wp, 0.023571_wp, 0.020714_wp, 0.021092_wp, 0.01951_wp, &
                      0.022279_wp, 0.022369_wp, 0.024528_wp, 0.023632_wp, 0.030801_wp, &
                      0.034516_wp, 0.033987_wp, 0.034206_wp, 0.032071_wp]
         case ("resp")
            charges = [ &
                      -0.389105_wp, 0.087886_wp, 0.155591_wp, -0.562255_wp, 0.268907_wp, &
                      -0.305056_wp, -0.071894_wp, -0.182209_wp, -0.113388_wp, -0.244008_wp, &
                      0.080274_wp, 0.098518_wp, 0.095375_wp, 0.011946_wp, -0.001353_wp, &
                      0.02223_wp, 0.123278_wp, 0.136527_wp, 0.13748_wp, 0.130193_wp, &
                      0.121776_wp, 0.1333_wp, 0.124882_wp, 0.141102_wp]
         case ("mbis")
            charges = [ &
                      -0.425062_wp, -0.181981_wp, -0.047464_wp, -0.453769_wp, 0.097261_wp, &
                      -0.202987_wp, -0.116931_wp, -0.153986_wp, -0.128606_wp, -0.178801_wp, &
                      0.125561_wp, 0.129422_wp, 0.145937_wp, 0.105467_wp, 0.10497_wp, &
                      0.106576_wp, 0.135414_wp, 0.143229_wp, 0.133436_wp, 0.131701_wp, &
                      0.131051_wp, 0.134984_wp, 0.13214_wp, 0.132452_wp]
         case ("chelpg")
            charges = [ &
                      -0.404754_wp, 0.098607_wp, 0.156236_wp, -0.579676_wp, 0.29715_wp, &
                      -0.331042_wp, -0.052473_wp, -0.200118_wp, -0.096131_wp, -0.266947_wp, &
                      0.083457_wp, 0.102261_wp, 0.098696_wp, 0.009976_wp, -0.003947_wp, &
                      0.02089_wp, 0.126978_wp, 0.140639_wp, 0.142022_wp, 0.135149_wp, &
                      0.118816_wp, 0.136555_wp, 0.122428_wp, 0.145228_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.075253_wp, -0.035091_wp, -0.003175_wp, -0.072295_wp, 0.009099_wp, &
                      -0.043529_wp, -0.041327_wp, -0.04681_wp, -0.040889_wp, -0.040615_wp, &
                      0.021977_wp, 0.024275_wp, 0.018242_wp, 0.022781_wp, 0.021091_wp, &
                      0.025013_wp, 0.02356_wp, 0.023793_wp, 0.025179_wp, 0.034586_wp, &
                      0.037979_wp, 0.037283_wp, 0.037713_wp, 0.036428_wp]
         case ("resp")
            charges = [ &
                      -0.40674_wp, 0.096283_wp, 0.166137_wp, -0.586535_wp, 0.270146_wp, &
                      -0.321453_wp, -0.081102_wp, -0.196708_wp, -0.122311_wp, -0.260073_wp, &
                      0.084913_wp, 0.103038_wp, 0.095177_wp, 0.012766_wp, -0.000397_wp, &
                      0.026716_wp, 0.129161_wp, 0.139429_wp, 0.145368_wp, 0.141399_wp, &
                      0.13182_wp, 0.143804_wp, 0.134644_wp, 0.154521_wp]
         case ("mbis")
            charges = [ &
                      -0.431414_wp, -0.184846_wp, -0.046167_wp, -0.463796_wp, 0.099022_wp, &
                      -0.21556_wp, -0.127522_wp, -0.165894_wp, -0.13938_wp, -0.190074_wp, &
                      0.128746_wp, 0.132182_wp, 0.142822_wp, 0.108897_wp, 0.108325_wp, &
                      0.111726_wp, 0.139311_wp, 0.143615_wp, 0.137938_wp, 0.141797_wp, &
                      0.140615_wp, 0.144497_wp, 0.141733_wp, 0.143441_wp]
         case ("chelpg")
            charges = [ &
                      -0.422671_wp, 0.107281_wp, 0.166698_wp, -0.604055_wp, 0.298153_wp, &
                      -0.34708_wp, -0.062099_wp, -0.21425_wp, -0.10538_wp, -0.282753_wp, &
                      0.088153_wp, 0.106844_wp, 0.098558_wp, 0.010744_wp, -0.003047_wp, &
                      0.0254_wp, 0.132884_wp, 0.143579_wp, 0.149938_wp, 0.14629_wp, &
                      0.128951_wp, 0.146992_wp, 0.132251_wp, 0.158617_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076486_wp, -0.03492_wp, -0.002934_wp, -0.073173_wp, 0.007389_wp, &
                      -0.046766_wp, -0.04555_wp, -0.051368_wp, -0.0449_wp, -0.043262_wp, &
                      0.022914_wp, 0.024919_wp, 0.015455_wp, 0.02465_wp, 0.022558_wp, &
                      0.027567_wp, 0.024636_wp, 0.022713_wp, 0.026847_wp, 0.037949_wp, &
                      0.040828_wp, 0.039934_wp, 0.040666_wp, 0.040347_wp]
         case ("resp")
            charges = [ &
                      -0.422135_wp, 0.105395_wp, 0.175273_wp, -0.60611_wp, 0.2694_wp, &
                      -0.334953_wp, -0.089288_wp, -0.209686_wp, -0.129628_wp, -0.273994_wp, &
                      0.089238_wp, 0.107015_wp, 0.093778_wp, 0.013778_wp, 0.000195_wp, &
                      0.030839_wp, 0.134022_wp, 0.140764_wp, 0.152531_wp, 0.151185_wp, &
                      0.140329_wp, 0.152678_wp, 0.142898_wp, 0.166477_wp]
         case ("mbis")
            charges = [ &
                      -0.436063_wp, -0.187103_wp, -0.045298_wp, -0.47154_wp, 0.100152_wp, &
                      -0.226384_wp, -0.136753_wp, -0.176494_wp, -0.148569_wp, -0.199748_wp, &
                      0.131611_wp, 0.134515_wp, 0.138775_wp, 0.112474_wp, 0.111293_wp, &
                      0.116516_wp, 0.142603_wp, 0.142945_wp, 0.142257_wp, 0.150673_wp, &
                      0.148665_wp, 0.152457_wp, 0.149863_wp, 0.153168_wp]
         case ("chelpg")
            charges = [ &
                      -0.438308_wp, 0.116661_wp, 0.175733_wp, -0.62369_wp, 0.29723_wp, &
                      -0.360292_wp, -0.070617_wp, -0.22694_wp, -0.112949_wp, -0.296478_wp, &
                      0.092527_wp, 0.110873_wp, 0.097211_wp, 0.011703_wp, -0.002511_wp, &
                      0.029549_wp, 0.137759_wp, 0.144943_wp, 0.15712_wp, 0.156024_wp, &
                      0.137532_wp, 0.155814_wp, 0.140553_wp, 0.170552_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (80)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071531_wp, 0.029678_wp, -0.072256_wp, -0.072236_wp, 0.008494_wp, &
                      -0.042834_wp, -0.037756_wp, -0.041821_wp, -0.037338_wp, -0.038293_wp, &
                      0.024576_wp, 0.024461_wp, 0.024462_wp, 0.022745_wp, 0.022537_wp, &
                      0.024117_wp, 0.024076_wp, 0.022585_wp, 0.022696_wp, 0.030303_wp, &
                      0.033704_wp, 0.034115_wp, 0.033952_wp, 0.031585_wp]
         case ("resp")
            charges = [ &
                      -0.495971_wp, 0.690975_wp, -0.521711_wp, -0.533735_wp, 0.081387_wp, &
                      -0.274324_wp, -0.082845_wp, -0.157017_wp, -0.100248_wp, -0.240092_wp, &
                      0.110792_wp, 0.102134_wp, 0.100743_wp, 0.115949_wp, 0.107988_wp, &
                      0.108071_wp, 0.110794_wp, 0.111769_wp, 0.119133_wp, 0.144693_wp, &
                      0.117715_wp, 0.125191_wp, 0.120346_wp, 0.138261_wp]
         case ("mbis")
            charges = [ &
                      -0.470666_wp, 0.159102_wp, -0.456447_wp, -0.456465_wp, 0.080598_wp, &
                      -0.199374_wp, -0.123648_wp, -0.149562_wp, -0.130476_wp, -0.174798_wp, &
                      0.141204_wp, 0.139455_wp, 0.139538_wp, 0.135188_wp, 0.136386_wp, &
                      0.142441_wp, 0.142438_wp, 0.136517_wp, 0.13511_wp, 0.145485_wp, &
                      0.130481_wp, 0.13423_wp, 0.131551_wp, 0.131734_wp]
         case ("chelpg")
            charges = [ &
                      -0.518395_wp, 0.703177_wp, -0.54196_wp, -0.553757_wp, 0.10265_wp, &
                      -0.297018_wp, -0.065148_wp, -0.172707_wp, -0.085308_wp, -0.256711_wp, &
                      0.116771_wp, 0.107605_wp, 0.106239_wp, 0.121324_wp, 0.113193_wp, &
                      0.112766_wp, 0.115456_wp, 0.116896_wp, 0.124422_wp, 0.149395_wp, &
                      0.114844_wp, 0.127966_wp, 0.118069_wp, 0.140234_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072189_wp, 0.029965_wp, -0.073464_wp, -0.073448_wp, 0.006209_wp, &
                      -0.04643_wp, -0.042541_wp, -0.046916_wp, -0.042061_wp, -0.041773_wp, &
                      0.026001_wp, 0.025487_wp, 0.02547_wp, 0.024508_wp, 0.023823_wp, &
                      0.023575_wp, 0.023517_wp, 0.023851_wp, 0.024465_wp, 0.0337_wp, &
                      0.037193_wp, 0.037444_wp, 0.037494_wp, 0.036141_wp]
         case ("resp")
            charges = [ &
                      -0.513276_wp, 0.725769_wp, -0.543612_wp, -0.555385_wp, 0.075237_wp, &
                      -0.286802_wp, -0.092486_wp, -0.171983_wp, -0.107748_wp, -0.25837_wp, &
                      0.115433_wp, 0.105973_wp, 0.104568_wp, 0.12303_wp, 0.113554_wp, &
                      0.109419_wp, 0.112088_wp, 0.117246_wp, 0.126133_wp, 0.155374_wp, &
                      0.127622_wp, 0.135798_wp, 0.129986_wp, 0.152433_wp]
         case ("mbis")
            charges = [ &
                      -0.479885_wp, 0.16795_wp, -0.467818_wp, -0.467821_wp, 0.079897_wp, &
                      -0.209516_wp, -0.13509_wp, -0.161412_wp, -0.141282_wp, -0.187045_wp, &
                      0.145105_wp, 0.142542_wp, 0.142591_wp, 0.140303_wp, 0.140583_wp, &
                      0.143275_wp, 0.143242_wp, 0.140684_wp, 0.14024_wp, 0.154681_wp, &
                      0.14021_wp, 0.14382_wp, 0.141259_wp, 0.143507_wp]
         case ("chelpg")
            charges = [ &
                      -0.535808_wp, 0.738226_wp, -0.564033_wp, -0.575567_wp, 0.095875_wp, &
                      -0.308916_wp, -0.075302_wp, -0.18722_wp, -0.093266_wp, -0.274448_wp, &
                      0.121432_wp, 0.111471_wp, 0.110091_wp, 0.128439_wp, 0.118805_wp, &
                      0.11416_wp, 0.116792_wp, 0.122416_wp, 0.131454_wp, 0.159964_wp, &
                      0.12485_wp, 0.138495_wp, 0.127789_wp, 0.1543_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072631_wp, 0.030161_wp, -0.074352_wp, -0.074341_wp, 0.004101_wp, &
                      -0.049553_wp, -0.046746_wp, -0.051467_wp, -0.046152_wp, -0.0447_wp, &
                      0.027481_wp, 0.02627_wp, 0.026232_wp, 0.02638_wp, 0.025018_wp, &
                      0.022708_wp, 0.022633_wp, 0.025028_wp, 0.026346_wp, 0.036653_wp, &
                      0.040104_wp, 0.040143_wp, 0.040496_wp, 0.040209_wp]
         case ("resp")
            charges = [ &
                      -0.526234_wp, 0.75395_wp, -0.560521_wp, -0.572107_wp, 0.068719_wp, &
                      -0.297172_wp, -0.100462_wp, -0.185687_wp, -0.113488_wp, -0.274358_wp, &
                      0.119416_wp, 0.108804_wp, 0.107374_wp, 0.129454_wp, 0.118227_wp, &
                      0.109338_wp, 0.111962_wp, 0.121844_wp, 0.132502_wp, 0.164537_wp, &
                      0.135957_wp, 0.144825_wp, 0.138095_wp, 0.165023_wp]
         case ("mbis")
            charges = [ &
                      -0.487009_wp, 0.174879_wp, -0.47673_wp, -0.476723_wp, 0.078837_wp, &
                      -0.21819_wp, -0.144866_wp, -0.171926_wp, -0.150427_wp, -0.197423_wp, &
                      0.148758_wp, 0.14486_wp, 0.14487_wp, 0.145185_wp, 0.144206_wp, &
                      0.143071_wp, 0.143008_wp, 0.144278_wp, 0.145139_wp, 0.162635_wp, &
                      0.148413_wp, 0.151856_wp, 0.149495_wp, 0.153826_wp]
         case ("chelpg")
            charges = [ &
                      -0.548843_wp, 0.766695_wp, -0.581083_wp, -0.592423_wp, 0.088698_wp, &
                      -0.318726_wp, -0.083736_wp, -0.200518_wp, -0.099423_wp, -0.289898_wp, &
                      0.125421_wp, 0.11432_wp, 0.112915_wp, 0.134886_wp, 0.123514_wp, &
                      0.114116_wp, 0.116702_wp, 0.127048_wp, 0.137843_wp, 0.169018_wp, &
                      0.133269_wp, 0.147454_wp, 0.135967_wp, 0.166782_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (81)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.066088_wp, -0.032978_wp, -0.033111_wp]
         case ("resp")
            charges = [ &
                      -0.103521_wp, 0.051881_wp, 0.051641_wp]
         case ("mbis")
            charges = [ &
                      0.038681_wp, -0.019254_wp, -0.019427_wp]
         case ("chelpg")
            charges = [ &
                      -0.103976_wp, 0.052108_wp, 0.051868_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.065137_wp, -0.032482_wp, -0.032655_wp]
         case ("resp")
            charges = [ &
                      -0.115848_wp, 0.058068_wp, 0.05778_wp]
         case ("mbis")
            charges = [ &
                      0.035262_wp, -0.01752_wp, -0.017742_wp]
         case ("chelpg")
            charges = [ &
                      -0.116331_wp, 0.058309_wp, 0.058022_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.064225_wp, -0.032007_wp, -0.032218_wp]
         case ("resp")
            charges = [ &
                      -0.125034_wp, 0.062685_wp, 0.062349_wp]
         case ("mbis")
            charges = [ &
                      0.032207_wp, -0.015968_wp, -0.016239_wp]
         case ("chelpg")
            charges = [ &
                      -0.125538_wp, 0.062937_wp, 0.062601_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (82)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.176415_wp, -0.04413_wp, -0.044099_wp, -0.04408_wp, -0.044108_wp]
         case ("resp")
            charges = [ &
                      -0.361901_wp, 0.090068_wp, 0.09085_wp, 0.090249_wp, 0.090734_wp]
         case ("mbis")
            charges = [ &
                      0.157619_wp, -0.039437_wp, -0.039398_wp, -0.039374_wp, -0.03941_wp]
         case ("chelpg")
            charges = [ &
                      -0.364848_wp, 0.090806_wp, 0.091588_wp, 0.090985_wp, 0.09147_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.175645_wp, -0.043947_wp, -0.043909_wp, -0.043878_wp, -0.043911_wp]
         case ("resp")
            charges = [ &
                      -0.37973_wp, 0.094488_wp, 0.095297_wp, 0.094755_wp, 0.095191_wp]
         case ("mbis")
            charges = [ &
                      0.155574_wp, -0.038938_wp, -0.038891_wp, -0.038851_wp, -0.038895_wp]
         case ("chelpg")
            charges = [ &
                      -0.382717_wp, 0.095235_wp, 0.096045_wp, 0.0955_wp, 0.095936_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.174894_wp, -0.04377_wp, -0.043726_wp, -0.043678_wp, -0.043721_wp]
         case ("resp")
            charges = [ &
                      -0.395891_wp, 0.098489_wp, 0.099324_wp, 0.098849_wp, 0.099229_wp]
         case ("mbis")
            charges = [ &
                      0.153424_wp, -0.038415_wp, -0.038358_wp, -0.038297_wp, -0.038355_wp]
         case ("chelpg")
            charges = [ &
                      -0.398934_wp, 0.09925_wp, 0.100086_wp, 0.099609_wp, 0.099989_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (83)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.034615_wp, -0.027055_wp, -0.035754_wp, 0.039579_wp, -0.035754_wp, &
                      -0.027055_wp, 0.039427_wp, 0.041434_wp, 0.045071_wp, 0.045071_wp, &
                      0.041434_wp, -0.091785_wp]
         case ("resp")
            charges = [ &
                      -0.112619_wp, -0.155592_wp, -0.073166_wp, 0.047255_wp, -0.075171_wp, &
                      -0.154489_wp, 0.127969_wp, 0.137647_wp, 0.123117_wp, 0.124146_wp, &
                      0.137234_wp, -0.126331_wp]
         case ("mbis")
            charges = [ &
                      -0.148235_wp, -0.118705_wp, -0.179051_wp, 0.161032_wp, -0.179047_wp, &
                      -0.118701_wp, 0.141813_wp, 0.140573_wp, 0.165791_wp, 0.16579_wp, &
                      0.140571_wp, -0.171832_wp]
         case ("chelpg")
            charges = [ &
                      -0.116719_wp, -0.151299_wp, -0.079732_wp, 0.054575_wp, -0.081739_wp, &
                      -0.150178_wp, 0.128699_wp, 0.137163_wp, 0.124486_wp, 0.125511_wp, &
                      0.136743_wp, -0.127509_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036154_wp, -0.029104_wp, -0.039407_wp, 0.031899_wp, -0.039407_wp, &
                      -0.029103_wp, 0.04712_wp, 0.049493_wp, 0.051095_wp, 0.051096_wp, &
                      0.049493_wp, -0.107022_wp]
         case ("resp")
            charges = [ &
                      -0.127387_wp, -0.163923_wp, -0.089955_wp, 0.047555_wp, -0.092053_wp, &
                      -0.162733_wp, 0.145639_wp, 0.155263_wp, 0.139806_wp, 0.140861_wp, &
                      0.154827_wp, -0.147899_wp]
         case ("mbis")
            charges = [ &
                      -0.158386_wp, -0.128261_wp, -0.191414_wp, 0.15409_wp, -0.191409_wp, &
                      -0.128257_wp, 0.157869_wp, 0.157227_wp, 0.179816_wp, 0.179815_wp, &
                      0.157225_wp, -0.188317_wp]
         case ("chelpg")
            charges = [ &
                      -0.132662_wp, -0.158412_wp, -0.097832_wp, 0.056085_wp, -0.099917_wp, &
                      -0.15721_wp, 0.146579_wp, 0.154562_wp, 0.141431_wp, 0.142478_wp, &
                      0.15412_wp, -0.14922_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036176_wp, -0.029369_wp, -0.040444_wp, 0.029737_wp, -0.040444_wp, &
                      -0.029368_wp, 0.049202_wp, 0.051601_wp, 0.05243_wp, 0.05243_wp, &
                      0.0516_wp, -0.111201_wp]
         case ("resp")
            charges = [ &
                      -0.130626_wp, -0.165283_wp, -0.094769_wp, 0.047659_wp, -0.09689_wp, &
                      -0.16407_wp, 0.150136_wp, 0.159655_wp, 0.143858_wp, 0.14492_wp, &
                      0.159213_wp, -0.153804_wp]
         case ("mbis")
            charges = [ &
                      -0.160346_wp, -0.130148_wp, -0.194688_wp, 0.151937_wp, -0.194683_wp, &
                      -0.130144_wp, 0.16194_wp, 0.161391_wp, 0.183065_wp, 0.183064_wp, &
                      0.161389_wp, -0.192778_wp]
         case ("chelpg")
            charges = [ &
                      -0.13616_wp, -0.159487_wp, -0.102958_wp, 0.056483_wp, -0.105064_wp, &
                      -0.158263_wp, 0.151125_wp, 0.158902_wp, 0.145545_wp, 0.146598_wp, &
                      0.158454_wp, -0.155175_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (84)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.123477_wp, 0.063511_wp, -0.06239_wp, -0.06231_wp, -0.062289_wp]
         case ("resp")
            charges = [ &
                      -0.299261_wp, 0.271011_wp, 0.009088_wp, 0.009469_wp, 0.009692_wp]
         case ("mbis")
            charges = [ &
                      0.064931_wp, 0.207026_wp, -0.090728_wp, -0.090636_wp, -0.090594_wp]
         case ("chelpg")
            charges = [ &
                      -0.302364_wp, 0.272103_wp, 0.009759_wp, 0.01014_wp, 0.010362_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.128017_wp, 0.0787_wp, -0.068984_wp, -0.068884_wp, -0.06885_wp]
         case ("resp")
            charges = [ &
                      -0.308236_wp, 0.311834_wp, -0.001526_wp, -0.001121_wp, -0.000951_wp]
         case ("mbis")
            charges = [ &
                      0.072877_wp, 0.230711_wp, -0.101293_wp, -0.101178_wp, -0.101119_wp]
         case ("chelpg")
            charges = [ &
                      -0.311196_wp, 0.312879_wp, -0.000887_wp, -0.000482_wp, -0.000314_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.12964_wp, 0.083655_wp, -0.071183_wp, -0.071076_wp, -0.071038_wp]
         case ("resp")
            charges = [ &
                      -0.311295_wp, 0.326112_wp, -0.005262_wp, -0.004851_wp, -0.004705_wp]
         case ("mbis")
            charges = [ &
                      0.075992_wp, 0.238261_wp, -0.104855_wp, -0.104732_wp, -0.104667_wp]
         case ("chelpg")
            charges = [ &
                      -0.314202_wp, 0.327142_wp, -0.004635_wp, -0.004225_wp, -0.00408_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (85)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.149498_wp, 0.023961_wp, -0.002804_wp, -0.030493_wp, -0.031337_wp, &
                      -0.030354_wp, -0.031391_wp, -0.030542_wp, 0.044634_wp, 0.044604_wp, &
                      0.038786_wp, 0.038662_wp, 0.038409_wp, 0.038646_wp, 0.038725_wp]
         case ("resp")
            charges = [ &
                      -0.148795_wp, -0.379551_wp, 0.301449_wp, -0.244062_wp, -0.105626_wp, &
                      -0.135373_wp, -0.10855_wp, -0.248618_wp, 0.179363_wp, 0.181442_wp, &
                      0.159458_wp, 0.128182_wp, 0.129371_wp, 0.130097_wp, 0.161214_wp]
         case ("mbis")
            charges = [ &
                      -0.224408_wp, -0.159317_wp, 0.083399_wp, -0.157876_wp, -0.131383_wp, &
                      -0.129534_wp, -0.131338_wp, -0.158159_wp, 0.157798_wp, 0.157741_wp, &
                      0.140959_wp, 0.137432_wp, 0.136676_wp, 0.137407_wp, 0.140611_wp]
         case ("chelpg")
            charges = [ &
                      -0.146637_wp, -0.397091_wp, 0.328453_wp, -0.264788_wp, -0.090121_wp, &
                      -0.150306_wp, -0.093355_wp, -0.268968_wp, 0.183178_wp, 0.18525_wp, &
                      0.163538_wp, 0.12583_wp, 0.132013_wp, 0.127838_wp, 0.165165_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.18029_wp, 0.029795_wp, -0.008083_wp, -0.035084_wp, -0.036603_wp, &
                      -0.034839_wp, -0.036653_wp, -0.035025_wp, 0.055049_wp, 0.055014_wp, &
                      0.045948_wp, 0.044981_wp, 0.044842_wp, 0.044969_wp, 0.045987_wp]
         case ("resp")
            charges = [ &
                      -0.204006_wp, -0.362841_wp, 0.295487_wp, -0.261875_wp, -0.122741_wp, &
                      -0.14827_wp, -0.126197_wp, -0.26482_wp, 0.199689_wp, 0.201603_wp, &
                      0.177283_wp, 0.145216_wp, 0.145735_wp, 0.147166_wp, 0.178571_wp]
         case ("mbis")
            charges = [ &
                      -0.265041_wp, -0.150238_wp, 0.07635_wp, -0.173538_wp, -0.146323_wp, &
                      -0.14187_wp, -0.146343_wp, -0.173652_wp, 0.174002_wp, 0.173938_wp, &
                      0.157662_wp, 0.152876_wp, 0.151876_wp, 0.152863_wp, 0.157446_wp]
         case ("chelpg")
            charges = [ &
                      -0.201847_wp, -0.380373_wp, 0.322163_wp, -0.282191_wp, -0.107719_wp, &
                      -0.162735_wp, -0.111476_wp, -0.284765_wp, 0.203513_wp, 0.20542_wp, &
                      0.181292_wp, 0.142966_wp, 0.148294_wp, 0.145006_wp, 0.182451_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.187191_wp, 0.031215_wp, -0.009135_wp, -0.036001_wp, -0.037664_wp, &
                      -0.035721_wp, -0.037709_wp, -0.035914_wp, 0.057349_wp, 0.057314_wp, &
                      0.047401_wp, 0.046237_wp, 0.046131_wp, 0.046229_wp, 0.047467_wp]
         case ("resp")
            charges = [ &
                      -0.217087_wp, -0.358815_wp, 0.293491_wp, -0.265046_wp, -0.126241_wp, &
                      -0.150929_wp, -0.129845_wp, -0.267492_wp, 0.204622_wp, 0.206477_wp, &
                      0.180687_wp, 0.148661_wp, 0.149071_wp, 0.150616_wp, 0.181831_wp]
         case ("mbis")
            charges = [ &
                      -0.274231_wp, -0.147683_wp, 0.074756_wp, -0.176635_wp, -0.149338_wp, &
                      -0.144264_wp, -0.149365_wp, -0.176706_wp, 0.177452_wp, 0.177388_wp, &
                      0.161016_wp, 0.155946_wp, 0.154897_wp, 0.155938_wp, 0.160837_wp]
         case ("chelpg")
            charges = [ &
                      -0.214927_wp, -0.376346_wp, 0.320105_wp, -0.285285_wp, -0.111307_wp, &
                      -0.16531_wp, -0.115211_wp, -0.28736_wp, 0.208448_wp, 0.210296_wp, &
                      0.184682_wp, 0.146429_wp, 0.151615_wp, 0.148474_wp, 0.185698_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (86)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06968_wp, 0.007472_wp, -0.035467_wp, -0.03603_wp, -0.033936_wp, &
                      -0.037652_wp, 0.032939_wp, -0.081068_wp, 0.032549_wp, 0.032538_wp, &
                      0.032088_wp, 0.035832_wp, 0.038084_wp, 0.039372_wp, 0.042959_wp]
         case ("resp")
            charges = [ &
                      -0.573111_wp, 0.410979_wp, -0.324724_wp, -0.094649_wp, -0.205898_wp, &
                      -0.069049_wp, -0.084357_wp, -0.110837_wp, 0.157617_wp, 0.157359_wp, &
                      0.153029_wp, 0.163023_wp, 0.137807_wp, 0.149461_wp, 0.13335_wp]
         case ("mbis")
            charges = [ &
                      -0.44873_wp, 0.101409_wp, -0.181996_wp, -0.133906_wp, -0.140864_wp, &
                      -0.163031_wp, 0.093605_wp, -0.16992_wp, 0.159708_wp, 0.159729_wp, &
                      0.142781_wp, 0.139114_wp, 0.138552_wp, 0.140852_wp, 0.162696_wp]
         case ("chelpg")
            charges = [ &
                      -0.599624_wp, 0.436281_wp, -0.341417_wp, -0.0813_wp, -0.220014_wp, &
                      -0.055457_wp, -0.100398_wp, -0.109035_wp, 0.16344_wp, 0.163183_wp, &
                      0.159045_wp, 0.165905_wp, 0.135837_wp, 0.152152_wp, 0.131401_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070607_wp, 0.005377_wp, -0.037817_wp, -0.038548_wp, -0.037217_wp, &
                      -0.041748_wp, 0.026776_wp, -0.093727_wp, 0.033659_wp, 0.033643_wp, &
                      0.037975_wp, 0.043418_wp, 0.044715_wp, 0.046065_wp, 0.048037_wp]
         case ("resp")
            charges = [ &
                      -0.593781_wp, 0.417894_wp, -0.340918_wp, -0.109063_wp, -0.216965_wp, &
                      -0.085363_wp, -0.087051_wp, -0.129603_wp, 0.163957_wp, 0.163688_wp, &
                      0.167393_wp, 0.181609_wp, 0.15428_wp, 0.165765_wp, 0.148159_wp]
         case ("mbis")
            charges = [ &
                      -0.46166_wp, 0.103674_wp, -0.195458_wp, -0.143999_wp, -0.152211_wp, &
                      -0.175367_wp, 0.08685_wp, -0.184003_wp, 0.163493_wp, 0.163507_wp, &
                      0.155268_wp, 0.155711_wp, 0.153122_wp, 0.15578_wp, 0.175293_wp]
         case ("chelpg")
            charges = [ &
                      -0.620381_wp, 0.442595_wp, -0.356851_wp, -0.096606_wp, -0.230148_wp, &
                      -0.072727_wp, -0.10236_wp, -0.127881_wp, 0.16982_wp, 0.169552_wp, &
                      0.17345_wp, 0.184361_wp, 0.152485_wp, 0.168282_wp, 0.146411_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070848_wp, 0.004591_wp, -0.038221_wp, -0.038966_wp, -0.038067_wp, &
                      -0.043266_wp, 0.024497_wp, -0.098277_wp, 0.033822_wp, 0.033804_wp, &
                      0.040119_wp, 0.046082_wp, 0.04702_wp, 0.048277_wp, 0.049434_wp]
         case ("resp")
            charges = [ &
                      -0.598902_wp, 0.418773_wp, -0.345076_wp, -0.113142_wp, -0.219903_wp, &
                      -0.0913_wp, -0.08783_wp, -0.136337_wp, 0.165379_wp, 0.165105_wp, &
                      0.172222_wp, 0.187695_wp, 0.159631_wp, 0.170962_wp, 0.152721_wp]
         case ("mbis")
            charges = [ &
                      -0.465131_wp, 0.103849_wp, -0.199063_wp, -0.14664_wp, -0.155407_wp, &
                      -0.179629_wp, 0.084252_wp, -0.18902_wp, 0.164261_wp, 0.164273_wp, &
                      0.159625_wp, 0.161208_wp, 0.15787_wp, 0.16054_wp, 0.179012_wp]
         case ("chelpg")
            charges = [ &
                      -0.62552_wp, 0.443296_wp, -0.360793_wp, -0.100938_wp, -0.232815_wp, &
                      -0.078947_wp, -0.102922_wp, -0.134639_wp, 0.171252_wp, 0.170979_wp, &
                      0.178289_wp, 0.190411_wp, 0.157885_wp, 0.173429_wp, 0.151032_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (87)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.057949_wp, -0.030745_wp, -0.060001_wp, 0.012074_wp, -0.071356_wp, &
                      0.083165_wp, -0.067003_wp, 0.039772_wp, 0.031523_wp, 0.02906_wp, &
                      0.029211_wp, 0.037523_wp, -0.204503_wp, 0.031213_wp, 0.030659_wp, &
                      0.167364_wp]
         case ("resp")
            charges = [ &
                      -0.328459_wp, -0.083547_wp, -0.393649_wp, 0.487579_wp, -0.590232_wp, &
                      0.494427_wp, -0.546135_wp, 0.187209_wp, 0.138877_wp, 0.170588_wp, &
                      0.199249_wp, 0.14424_wp, -0.577278_wp, 0.146674_wp, 0.146985_wp, &
                      0.403474_wp]
         case ("mbis")
            charges = [ &
                      -0.31539_wp, -0.061975_wp, -0.277196_wp, 0.229414_wp, -0.427665_wp, &
                      0.428456_wp, -0.470532_wp, 0.170402_wp, 0.147664_wp, 0.144224_wp, &
                      0.155728_wp, 0.131128_wp, -0.586497_wp, 0.147295_wp, 0.146753_wp, &
                      0.438198_wp]
         case ("chelpg")
            charges = [ &
                      -0.359023_wp, -0.056851_wp, -0.424296_wp, 0.529758_wp, -0.630273_wp, &
                      0.528567_wp, -0.575249_wp, 0.193121_wp, 0.144676_wp, 0.175838_wp, &
                      0.207076_wp, 0.140216_wp, -0.583648_wp, 0.152728_wp, 0.153068_wp, &
                      0.404293_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071194_wp, -0.040071_wp, -0.06838_wp, 0.008022_wp, -0.070906_wp, &
                      0.074292_wp, -0.067943_wp, 0.041007_wp, 0.035047_wp, 0.034987_wp, &
                      0.04458_wp, 0.043437_wp, -0.220563_wp, 0.035327_wp, 0.035295_wp, &
                      0.187072_wp]
         case ("resp")
            charges = [ &
                      -0.373033_wp, -0.09428_wp, -0.429302_wp, 0.488502_wp, -0.608215_wp, &
                      0.51471_wp, -0.571786_wp, 0.202571_wp, 0.150652_wp, 0.191786_wp, &
                      0.232172_wp, 0.162206_wp, -0.639985_wp, 0.159143_wp, 0.159191_wp, &
                      0.455668_wp]
         case ("mbis")
            charges = [ &
                      -0.350399_wp, -0.07591_wp, -0.304288_wp, 0.231365_wp, -0.441887_wp, &
                      0.432762_wp, -0.488969_wp, 0.182502_wp, 0.156598_wp, 0.161908_wp, &
                      0.185973_wp, 0.147376_wp, -0.635803_wp, 0.156949_wp, 0.15731_wp, &
                      0.484522_wp]
         case ("chelpg")
            charges = [ &
                      -0.4034_wp, -0.067881_wp, -0.459736_wp, 0.530667_wp, -0.648217_wp, &
                      0.548766_wp, -0.601137_wp, 0.208446_wp, 0.156508_wp, 0.196994_wp, &
                      0.240003_wp, 0.158257_wp, -0.646362_wp, 0.165258_wp, 0.165336_wp, &
                      0.456497_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072779_wp, -0.041155_wp, -0.069238_wp, 0.007705_wp, -0.070637_wp, &
                      0.073352_wp, -0.067975_wp, 0.040915_wp, 0.035377_wp, 0.035577_wp, &
                      0.046353_wp, 0.043943_wp, -0.222169_wp, 0.035812_wp, 0.035851_wp, &
                      0.189075_wp]
         case ("resp")
            charges = [ &
                      -0.377866_wp, -0.095561_wp, -0.432714_wp, 0.488026_wp, -0.608953_wp, &
                      0.516341_wp, -0.573932_wp, 0.203861_wp, 0.151693_wp, 0.193897_wp, &
                      0.235633_wp, 0.163944_wp, -0.646234_wp, 0.160444_wp, 0.160489_wp, &
                      0.460934_wp]
         case ("mbis")
            charges = [ &
                      -0.354381_wp, -0.077433_wp, -0.306987_wp, 0.23151_wp, -0.442736_wp, &
                      0.433057_wp, -0.4907_wp, 0.183489_wp, 0.157408_wp, 0.163675_wp, &
                      0.18925_wp, 0.14892_wp, -0.640641_wp, 0.157996_wp, 0.158483_wp, &
                      0.4891_wp]
         case ("chelpg")
            charges = [ &
                      -0.40821_wp, -0.069204_wp, -0.463119_wp, 0.530181_wp, -0.648941_wp, &
                      0.550384_wp, -0.603296_wp, 0.209732_wp, 0.157553_wp, 0.199099_wp, &
                      0.243462_wp, 0.160004_wp, -0.65261_wp, 0.166564_wp, 0.166639_wp, &
                      0.461762_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (88)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.053967_wp, -0.037722_wp, -0.007091_wp, 0.078229_wp, -0.066018_wp, &
                      -0.038621_wp, 0.032649_wp, 0.033548_wp, -0.068874_wp, -0.198114_wp, &
                      0.03127_wp, 0.03581_wp, 0.168139_wp, 0.02881_wp, 0.030962_wp, &
                      0.030984_wp]
         case ("resp")
            charges = [ &
                      -0.192683_wp, -0.286454_wp, 0.285129_wp, 0.229465_wp, -0.314661_wp, &
                      -0.114394_wp, 0.147012_wp, 0.16918_wp, -0.532129_wp, -0.541615_wp, &
                      0.159088_wp, 0.14083_wp, 0.403828_wp, 0.131588_wp, 0.157788_wp, &
                      0.158027_wp]
         case ("mbis")
            charges = [ &
                      -0.197948_wp, -0.146335_wp, -0.001343_wp, 0.336206_wp, -0.34412_wp, &
                      -0.074713_wp, 0.141879_wp, 0.133772_wp, -0.420931_wp, -0.571885_wp, &
                      0.148072_wp, 0.130542_wp, 0.438343_wp, 0.136034_wp, 0.146308_wp, &
                      0.146116_wp]
         case ("chelpg")
            charges = [ &
                      -0.194308_wp, -0.289261_wp, 0.294722_wp, 0.230984_wp, -0.319662_wp, &
                      -0.111881_wp, 0.147654_wp, 0.169968_wp, -0.553691_wp, -0.54313_wp, &
                      0.160507_wp, 0.140574_wp, 0.404595_wp, 0.136776_wp, 0.162968_wp, &
                      0.163183_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.062881_wp, -0.045657_wp, -0.014238_wp, 0.072609_wp, -0.066775_wp, &
                      -0.045805_wp, 0.037155_wp, 0.039472_wp, -0.072039_wp, -0.209612_wp, &
                      0.044767_wp, 0.041792_wp, 0.185881_wp, 0.034352_wp, 0.030434_wp, &
                      0.03054_wp]
         case ("resp")
            charges = [ &
                      -0.215586_wp, -0.308769_wp, 0.268865_wp, 0.256359_wp, -0.338626_wp, &
                      -0.136011_wp, 0.163296_wp, 0.187886_wp, -0.543496_wp, -0.597134_wp, &
                      0.188925_wp, 0.158629_wp, 0.450454_wp, 0.146318_wp, 0.159282_wp, &
                      0.159608_wp]
         case ("mbis")
            charges = [ &
                      -0.218086_wp, -0.166283_wp, -0.009281_wp, 0.338347_wp, -0.355791_wp, &
                      -0.094768_wp, 0.156209_wp, 0.15017_wp, -0.435916_wp, -0.612228_wp, &
                      0.175037_wp, 0.146747_wp, 0.479179_wp, 0.149366_wp, 0.148671_wp, &
                      0.148622_wp]
         case ("chelpg")
            charges = [ &
                      -0.217092_wp, -0.311645_wp, 0.278514_wp, 0.257868_wp, -0.343542_wp, &
                      -0.133679_wp, 0.163919_wp, 0.188682_wp, -0.56515_wp, -0.598653_wp, &
                      0.190333_wp, 0.158429_wp, 0.451217_wp, 0.151532_wp, 0.164482_wp, &
                      0.164784_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.064759_wp, -0.047418_wp, -0.015904_wp, 0.071457_wp, -0.066578_wp, &
                      -0.047107_wp, 0.03797_wp, 0.040579_wp, -0.072707_wp, -0.211826_wp, &
                      0.047783_wp, 0.043047_wp, 0.189433_wp, 0.035541_wp, 0.03018_wp, &
                      0.030305_wp]
         case ("resp")
            charges = [ &
                      -0.220401_wp, -0.313368_wp, 0.264626_wp, 0.261537_wp, -0.342314_wp, &
                      -0.140285_wp, 0.166484_wp, 0.191566_wp, -0.544833_wp, -0.608103_wp, &
                      0.195129_wp, 0.162237_wp, 0.459737_wp, 0.149338_wp, 0.15915_wp, &
                      0.159499_wp]
         case ("mbis")
            charges = [ &
                      -0.222289_wp, -0.170442_wp, -0.011444_wp, 0.338519_wp, -0.357236_wp, &
                      -0.098692_wp, 0.158975_wp, 0.153381_wp, -0.438523_wp, -0.619964_wp, &
                      0.18073_wp, 0.15002_wp, 0.487126_wp, 0.152157_wp, 0.148847_wp, &
                      0.148832_wp]
         case ("chelpg")
            charges = [ &
                      -0.221892_wp, -0.316255_wp, 0.274277_wp, 0.263051_wp, -0.347216_wp, &
                      -0.137982_wp, 0.167105_wp, 0.192364_wp, -0.566488_wp, -0.609624_wp, &
                      0.196537_wp, 0.162046_wp, 0.460499_wp, 0.154555_wp, 0.164349_wp, &
                      0.164674_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (89)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.0385_wp, -0.03852_wp, -0.038495_wp, -0.038494_wp, -0.038519_wp, &
                      -0.038499_wp, 0.017616_wp, 0.017609_wp, 0.017619_wp, 0.017612_wp, &
                      0.017608_wp, 0.017622_wp, 0.02089_wp, 0.020894_wp, 0.020888_wp, &
                      0.020889_wp, 0.020894_wp, 0.020889_wp]
         case ("resp")
            charges = [ &
                      -0.00609_wp, -0.02502_wp, -0.01024_wp, -0.025012_wp, -0.017199_wp, &
                      -0.02582_wp, 0.014392_wp, 0.01805_wp, 0.015345_wp, 0.018235_wp, &
                      0.016413_wp, 0.018134_wp, -0.00178_wp, 0.001779_wp, -7.1e-05_wp, &
                      0.003034_wp, 0.002588_wp, 0.003263_wp]
         case ("mbis")
            charges = [ &
                      -0.204456_wp, -0.204453_wp, -0.204509_wp, -0.204457_wp, -0.204456_wp, &
                      -0.204511_wp, 0.101423_wp, 0.101444_wp, 0.101422_wp, 0.101401_wp, &
                      0.101441_wp, 0.101441_wp, 0.103035_wp, 0.10302_wp, 0.103082_wp, &
                      0.103034_wp, 0.10302_wp, 0.103083_wp]
         case ("chelpg")
            charges = [ &
                      -0.004703_wp, -0.027591_wp, -0.009167_wp, -0.027375_wp, -0.016661_wp, &
                      -0.028305_wp, 0.0144_wp, 0.018799_wp, 0.015428_wp, 0.018982_wp, &
                      0.016637_wp, 0.018896_wp, -0.001862_wp, 0.002358_wp, -5.5e-05_wp, &
                      0.003598_wp, 0.002774_wp, 0.003848_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.039081_wp, -0.039105_wp, -0.039081_wp, -0.039076_wp, -0.039104_wp, &
                      -0.039086_wp, 0.018024_wp, 0.018017_wp, 0.018021_wp, 0.01802_wp, &
                      0.018016_wp, 0.018025_wp, 0.02107_wp, 0.021071_wp, 0.021067_wp, &
                      0.021069_wp, 0.021071_wp, 0.021067_wp]
         case ("resp")
            charges = [ &
                      -0.006022_wp, -0.025275_wp, -0.010321_wp, -0.025273_wp, -0.017433_wp, &
                      -0.026209_wp, 0.015143_wp, 0.018861_wp, 0.016105_wp, 0.019051_wp, &
                      0.017235_wp, 0.018966_wp, -0.002451_wp, 0.001155_wp, -0.000701_wp, &
                      0.002452_wp, 0.002014_wp, 0.002703_wp]
         case ("mbis")
            charges = [ &
                      -0.206468_wp, -0.206486_wp, -0.20654_wp, -0.206467_wp, -0.206489_wp, &
                      -0.206544_wp, 0.102661_wp, 0.102687_wp, 0.102653_wp, 0.102638_wp, &
                      0.102684_wp, 0.102674_wp, 0.103822_wp, 0.103806_wp, 0.103871_wp, &
                      0.103821_wp, 0.103806_wp, 0.103871_wp]
         case ("chelpg")
            charges = [ &
                      -0.004607_wp, -0.027878_wp, -0.009232_wp, -0.027666_wp, -0.01689_wp, &
                      -0.028734_wp, 0.015149_wp, 0.01962_wp, 0.016188_wp, 0.019809_wp, &
                      0.017463_wp, 0.01974_wp, -0.002536_wp, 0.00174_wp, -0.000685_wp, &
                      0.003022_wp, 0.002203_wp, 0.003297_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.039663_wp, -0.039692_wp, -0.039668_wp, -0.039657_wp, -0.039691_wp, &
                      -0.039674_wp, 0.018488_wp, 0.018482_wp, 0.018479_wp, 0.018483_wp, &
                      0.018481_wp, 0.018483_wp, 0.021195_wp, 0.021193_wp, 0.02119_wp, &
                      0.021194_wp, 0.021192_wp, 0.02119_wp]
         case ("resp")
            charges = [ &
                      -0.005711_wp, -0.025317_wp, -0.010176_wp, -0.025314_wp, -0.017452_wp, &
                      -0.026396_wp, 0.015923_wp, 0.019708_wp, 0.016895_wp, 0.019901_wp, &
                      0.018094_wp, 0.019835_wp, -0.003379_wp, 0.000275_wp, -0.001586_wp, &
                      0.001617_wp, 0.001188_wp, 0.001893_wp]
         case ("mbis")
            charges = [ &
                      -0.208486_wp, -0.208528_wp, -0.208578_wp, -0.208482_wp, -0.20853_wp, &
                      -0.208584_wp, 0.10401_wp, 0.104041_wp, 0.103994_wp, 0.103985_wp, &
                      0.104038_wp, 0.104017_wp, 0.104507_wp, 0.104489_wp, 0.104558_wp, &
                      0.104507_wp, 0.104489_wp, 0.104557_wp]
         case ("chelpg")
            charges = [ &
                      -0.004255_wp, -0.02795_wp, -0.009061_wp, -0.027731_wp, -0.016893_wp, &
                      -0.028959_wp, 0.015922_wp, 0.020473_wp, 0.016975_wp, 0.020665_wp, &
                      0.018322_wp, 0.020617_wp, -0.003474_wp, 0.000865_wp, -0.001575_wp, &
                      0.00219_wp, 0.001376_wp, 0.002494_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (90)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.172867_wp, -0.04277_wp, -0.032514_wp, -0.037271_wp, -0.032516_wp, &
                      -0.042775_wp, -0.28323_wp, 0.034796_wp, 0.025239_wp, 0.021383_wp, &
                      0.025241_wp, 0.034789_wp, 0.036632_wp, 0.028115_wp, 0.027267_wp, &
                      0.028118_wp, 0.036629_wp]
         case ("resp")
            charges = [ &
                      0.499212_wp, -0.145169_wp, 0.028471_wp, -0.102879_wp, 0.028612_wp, &
                      -0.15178_wp, -0.513677_wp, 0.051499_wp, 0.024342_wp, 0.04003_wp, &
                      0.025008_wp, 0.053632_wp, 0.055108_wp, 0.012515_wp, 0.02506_wp, &
                      0.012925_wp, 0.057089_wp]
         case ("mbis")
            charges = [ &
                      0.623756_wp, -0.406522_wp, -0.169135_wp, -0.238093_wp, -0.169152_wp, &
                      -0.406531_wp, -0.522741_wp, 0.147278_wp, 0.110737_wp, 0.111955_wp, &
                      0.11075_wp, 0.147275_wp, 0.159515_wp, 0.111718_wp, 0.117955_wp, &
                      0.111723_wp, 0.159514_wp]
         case ("chelpg")
            charges = [ &
                      0.517035_wp, -0.169464_wp, 0.045819_wp, -0.121573_wp, 0.045926_wp, &
                      -0.176133_wp, -0.516781_wp, 0.056544_wp, 0.021587_wp, 0.044326_wp, &
                      0.022278_wp, 0.058682_wp, 0.060061_wp, 0.010642_wp, 0.027918_wp, &
                      0.011067_wp, 0.062067_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.180536_wp, -0.041714_wp, -0.030197_wp, -0.03586_wp, -0.0302_wp, &
                      -0.041721_wp, -0.336969_wp, 0.045952_wp, 0.027328_wp, 0.025869_wp, &
                      0.027327_wp, 0.045947_wp, 0.034391_wp, 0.032487_wp, 0.029947_wp, &
                      0.03249_wp, 0.034387_wp]
         case ("resp")
            charges = [ &
                      0.538675_wp, -0.147797_wp, 0.053852_wp, -0.118729_wp, 0.054082_wp, &
                      -0.154899_wp, -0.604518_wp, 0.068195_wp, 0.019921_wp, 0.04759_wp, &
                      0.020596_wp, 0.070465_wp, 0.048425_wp, 0.012731_wp, 0.027689_wp, &
                      0.013163_wp, 0.050559_wp]
         case ("mbis")
            charges = [ &
                      0.677542_wp, -0.42854_wp, -0.165344_wp, -0.245522_wp, -0.165367_wp, &
                      -0.428554_wp, -0.607232_wp, 0.169725_wp, 0.113326_wp, 0.120237_wp, &
                      0.113336_wp, 0.169727_wp, 0.158423_wp, 0.11851_wp, 0.122796_wp, &
                      0.118517_wp, 0.15842_wp]
         case ("chelpg")
            charges = [ &
                      0.557858_wp, -0.174357_wp, 0.074057_wp, -0.140047_wp, 0.07424_wp, &
                      -0.181506_wp, -0.607783_wp, 0.073672_wp, 0.016559_wp, 0.052326_wp, &
                      0.017263_wp, 0.075942_wp, 0.053751_wp, 0.010363_wp, 0.030941_wp, &
                      0.010814_wp, 0.055907_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.181084_wp, -0.041627_wp, -0.029969_wp, -0.035673_wp, -0.029972_wp, &
                      -0.041634_wp, -0.341591_wp, 0.047001_wp, 0.027418_wp, 0.026365_wp, &
                      0.027417_wp, 0.046996_wp, 0.034072_wp, 0.03292_wp, 0.030203_wp, &
                      0.032923_wp, 0.034067_wp]
         case ("resp")
            charges = [ &
                      0.54159_wp, -0.147941_wp, 0.056041_wp, -0.119843_wp, 0.056273_wp, &
                      -0.15508_wp, -0.61226_wp, 0.069842_wp, 0.019348_wp, 0.048372_wp, &
                      0.020024_wp, 0.072123_wp, 0.047628_wp, 0.012881_wp, 0.027913_wp, &
                      0.013315_wp, 0.049774_wp]
         case ("mbis")
            charges = [ &
                      0.681895_wp, -0.430342_wp, -0.164928_wp, -0.246078_wp, -0.164951_wp, &
                      -0.430356_wp, -0.61446_wp, 0.171802_wp, 0.113367_wp, 0.121112_wp, &
                      0.113376_wp, 0.171805_wp, 0.158091_wp, 0.119178_wp, 0.123216_wp, &
                      0.119185_wp, 0.158089_wp]
         case ("chelpg")
            charges = [ &
                      0.560871_wp, -0.17466_wp, 0.07646_wp, -0.141351_wp, 0.076645_wp, &
                      -0.181844_wp, -0.615537_wp, 0.075349_wp, 0.015938_wp, 0.05314_wp, &
                      0.016643_wp, 0.077629_wp, 0.052978_wp, 0.010474_wp, 0.031192_wp, &
                      0.010928_wp, 0.055145_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (91)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.038993_wp, -0.039617_wp, -0.040351_wp, -0.040124_wp, -0.039261_wp, &
                      0.019421_wp, 0.020975_wp, 0.021805_wp, 0.021566_wp, 0.018438_wp, &
                      0.018827_wp, 0.018516_wp, 0.019424_wp, 0.019015_wp, 0.020358_wp]
         case ("resp")
            charges = [ &
                      -0.02349_wp, -0.017545_wp, -0.01841_wp, -0.007994_wp, -0.011125_wp, &
                      0.006654_wp, -0.000924_wp, -0.005296_wp, -0.007605_wp, 0.012259_wp, &
                      0.01018_wp, 0.018825_wp, 0.02532_wp, 0.020599_wp, -0.001448_wp]
         case ("mbis")
            charges = [ &
                      -0.192129_wp, -0.205374_wp, -0.212791_wp, -0.211246_wp, -0.198933_wp, &
                      0.097142_wp, 0.102421_wp, 0.105722_wp, 0.104877_wp, 0.098733_wp, &
                      0.096705_wp, 0.101568_wp, 0.107721_wp, 0.105601_wp, 0.099981_wp]
         case ("chelpg")
            charges = [ &
                      -0.024528_wp, -0.017763_wp, -0.019261_wp, -0.008015_wp, -0.011017_wp, &
                      0.006997_wp, -0.000748_wp, -0.005041_wp, -0.007585_wp, 0.012341_wp, &
                      0.010537_wp, 0.019055_wp, 0.025694_wp, 0.020755_wp, -0.001422_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.039726_wp, -0.040321_wp, -0.041051_wp, -0.040823_wp, -0.03998_wp, &
                      0.019698_wp, 0.021212_wp, 0.022068_wp, 0.021817_wp, 0.018885_wp, &
                      0.019176_wp, 0.019003_wp, 0.019928_wp, 0.019519_wp, 0.020593_wp]
         case ("resp")
            charges = [ &
                      -0.023921_wp, -0.016516_wp, -0.017301_wp, -0.006509_wp, -0.010776_wp, &
                      0.006368_wp, -0.002073_wp, -0.006746_wp, -0.009076_wp, 0.012493_wp, &
                      0.010192_wp, 0.019155_wp, 0.025931_wp, 0.021048_wp, -0.002268_wp]
         case ("mbis")
            charges = [ &
                      -0.19454_wp, -0.207831_wp, -0.215297_wp, -0.213721_wp, -0.201386_wp, &
                      0.09822_wp, 0.103422_wp, 0.106773_wp, 0.105902_wp, 0.100139_wp, &
                      0.09792_wp, 0.103053_wp, 0.109243_wp, 0.107119_wp, 0.100983_wp]
         case ("chelpg")
            charges = [ &
                      -0.025032_wp, -0.016673_wp, -0.018132_wp, -0.00646_wp, -0.010658_wp, &
                      0.006728_wp, -0.001914_wp, -0.006508_wp, -0.009081_wp, 0.01257_wp, &
                      0.010566_wp, 0.019368_wp, 0.02629_wp, 0.02118_wp, -0.002246_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.040511_wp, -0.041071_wp, -0.041794_wp, -0.041565_wp, -0.040748_wp, &
                      0.019969_wp, 0.021422_wp, 0.022311_wp, 0.022044_wp, 0.019394_wp, &
                      0.019549_wp, 0.01957_wp, 0.020516_wp, 0.020108_wp, 0.020802_wp]
         case ("resp")
            charges = [ &
                      -0.024166_wp, -0.015032_wp, -0.015693_wp, -0.004472_wp, -0.010106_wp, &
                      0.005871_wp, -0.003604_wp, -0.008627_wp, -0.01098_wp, 0.012669_wp, &
                      0.010074_wp, 0.01945_wp, 0.026548_wp, 0.021479_wp, -0.00341_wp]
         case ("mbis")
            charges = [ &
                      -0.197102_wp, -0.210462_wp, -0.217966_wp, -0.216362_wp, -0.204006_wp, &
                      0.099324_wp, 0.104408_wp, 0.107821_wp, 0.106917_wp, 0.101707_wp, &
                      0.09922_wp, 0.104728_wp, 0.110962_wp, 0.108836_wp, 0.101972_wp]
         case ("chelpg")
            charges = [ &
                      -0.02536_wp, -0.015111_wp, -0.016492_wp, -0.004331_wp, -0.009972_wp, &
                      0.006249_wp, -0.003467_wp, -0.00841_wp, -0.011018_wp, 0.012738_wp, &
                      0.010466_wp, 0.019638_wp, 0.026886_wp, 0.021579_wp, -0.003395_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (92)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.242959_wp, 0.06378_wp, -0.047402_wp, -0.039777_wp, -0.037931_wp, &
                      -0.041682_wp, 0.146233_wp, 0.017802_wp, 0.021478_wp, 0.024174_wp, &
                      0.020611_wp, 0.020048_wp, 0.0236_wp, 0.019497_wp, 0.023414_wp, &
                      0.029114_wp]
         case ("resp")
            charges = [ &
                      -0.626835_wp, 0.342343_wp, -0.128911_wp, 0.02963_wp, -0.060741_wp, &
                      -0.02918_wp, 0.37516_wp, -0.054789_wp, 0.020201_wp, 0.035627_wp, &
                      0.006149_wp, 0.002061_wp, 0.036771_wp, 0.012033_wp, 0.007515_wp, &
                      0.032966_wp]
         case ("mbis")
            charges = [ &
                      -0.61207_wp, 0.232642_wp, -0.327285_wp, -0.178173_wp, -0.185819_wp, &
                      -0.277913_wp, 0.411599_wp, 0.047818_wp, 0.113656_wp, 0.127054_wp, &
                      0.103785_wp, 0.097158_wp, 0.107974_wp, 0.095375_wp, 0.113273_wp, &
                      0.130927_wp]
         case ("chelpg")
            charges = [ &
                      -0.629066_wp, 0.351278_wp, -0.138616_wp, 0.036965_wp, -0.067093_wp, &
                      -0.029785_wp, 0.375616_wp, -0.057253_wp, 0.02207_wp, 0.037711_wp, &
                      0.004861_wp, 0.000775_wp, 0.038289_wp, 0.01362_wp, 0.00753_wp, &
                      0.033097_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.267088_wp, 0.064375_wp, -0.045782_wp, -0.041954_wp, -0.041092_wp, &
                      -0.045112_wp, 0.161034_wp, 0.023395_wp, 0.028286_wp, 0.028851_wp, &
                      0.018194_wp, 0.022372_wp, 0.01917_wp, 0.02203_wp, 0.028878_wp, &
                      0.024441_wp]
         case ("resp")
            charges = [ &
                      -0.718831_wp, 0.384195_wp, -0.115368_wp, 0.006148_wp, -0.038499_wp, &
                      -0.059682_wp, 0.426963_wp, -0.052208_wp, 0.026963_wp, 0.037013_wp, &
                      0.005883_wp, 0.008969_wp, 0.023573_wp, 0.013543_wp, 0.022317_wp, &
                      0.029021_wp]
         case ("mbis")
            charges = [ &
                      -0.680075_wp, 0.250836_wp, -0.334553_wp, -0.184698_wp, -0.189538_wp, &
                      -0.294749_wp, 0.455311_wp, 0.05591_wp, 0.126195_wp, 0.135065_wp, &
                      0.100994_wp, 0.102872_wp, 0.101863_wp, 0.101763_wp, 0.126405_wp, &
                      0.126399_wp]
         case ("chelpg")
            charges = [ &
                      -0.720819_wp, 0.393593_wp, -0.123764_wp, 0.010421_wp, -0.041395_wp, &
                      -0.063401_wp, 0.427219_wp, -0.054722_wp, 0.028635_wp, 0.038852_wp, &
                      0.005232_wp, 0.008369_wp, 0.024402_wp, 0.014448_wp, 0.023084_wp, &
                      0.029844_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.268901_wp, 0.064442_wp, -0.04559_wp, -0.042099_wp, -0.041348_wp, &
                      -0.045379_wp, 0.162089_wp, 0.023874_wp, 0.028858_wp, 0.029282_wp, &
                      0.017942_wp, 0.022586_wp, 0.018687_wp, 0.022249_wp, 0.029365_wp, &
                      0.023942_wp]
         case ("resp")
            charges = [ &
                      -0.725794_wp, 0.387271_wp, -0.113957_wp, 0.004369_wp, -0.036851_wp, &
                      -0.062144_wp, 0.430779_wp, -0.051862_wp, 0.027489_wp, 0.037176_wp, &
                      0.005724_wp, 0.009563_wp, 0.022312_wp, 0.013767_wp, 0.023644_wp, &
                      0.028513_wp]
         case ("mbis")
            charges = [ &
                      -0.685121_wp, 0.252176_wp, -0.334947_wp, -0.185132_wp, -0.189756_wp, &
                      -0.295978_wp, 0.458439_wp, 0.056618_wp, 0.1272_wp, 0.135757_wp, &
                      0.100638_wp, 0.103359_wp, 0.101126_wp, 0.1023_wp, 0.127533_wp, &
                      0.125787_wp]
         case ("chelpg")
            charges = [ &
                      -0.727762_wp, 0.396696_wp, -0.122243_wp, 0.008397_wp, -0.039474_wp, &
                      -0.066087_wp, 0.43102_wp, -0.054378_wp, 0.029145_wp, 0.038993_wp, &
                      0.005124_wp, 0.009018_wp, 0.023087_wp, 0.014617_wp, 0.024464_wp, &
                      0.029385_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (93)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.28564_wp, 0.175036_wp, -0.044271_wp, -0.044269_wp, -0.036596_wp, &
                      -0.036581_wp, 0.041248_wp, 0.041406_wp, 0.041237_wp, 0.041392_wp, &
                      0.028431_wp, 0.025073_wp, 0.028437_wp, 0.025096_wp]
         case ("resp")
            charges = [ &
                      -0.537665_wp, 0.572902_wp, -0.211003_wp, -0.206714_wp, -0.046921_wp, &
                      -0.047642_wp, 0.087955_wp, 0.076914_wp, 0.086082_wp, 0.076264_wp, &
                      0.031486_wp, 0.043736_wp, 0.03167_wp, 0.042936_wp]
         case ("mbis")
            charges = [ &
                      -0.519536_wp, 0.620269_wp, -0.408405_wp, -0.408357_wp, -0.190699_wp, &
                      -0.190789_wp, 0.159383_wp, 0.162272_wp, 0.15934_wp, 0.16229_wp, &
                      0.114436_wp, 0.112628_wp, 0.114462_wp, 0.112705_wp]
         case ("chelpg")
            charges = [ &
                      -0.539352_wp, 0.579861_wp, -0.21951_wp, -0.215078_wp, -0.047135_wp, &
                      -0.047732_wp, 0.090195_wp, 0.079229_wp, 0.088269_wp, 0.078532_wp, &
                      0.032152_wp, 0.044556_wp, 0.03229_wp, 0.043722_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.338775_wp, 0.18238_wp, -0.042656_wp, -0.04266_wp, -0.033742_wp, &
                      -0.033731_wp, 0.048705_wp, 0.04406_wp, 0.048752_wp, 0.043997_wp, &
                      0.032638_wp, 0.029187_wp, 0.032649_wp, 0.029194_wp]
         case ("resp")
            charges = [ &
                      -0.630138_wp, 0.614358_wp, -0.204866_wp, -0.201446_wp, -0.039079_wp, &
                      -0.039737_wp, 0.095956_wp, 0.076593_wp, 0.094367_wp, 0.076164_wp, &
                      0.033654_wp, 0.045538_wp, 0.033967_wp, 0.04467_wp]
         case ("mbis")
            charges = [ &
                      -0.602668_wp, 0.672836_wp, -0.429719_wp, -0.429695_wp, -0.188734_wp, &
                      -0.188831_wp, 0.174727_wp, 0.169317_wp, 0.174788_wp, 0.169253_wp, &
                      0.120574_wp, 0.118744_wp, 0.120611_wp, 0.118797_wp]
         case ("chelpg")
            charges = [ &
                      -0.631819_wp, 0.621356_wp, -0.213526_wp, -0.209975_wp, -0.0391_wp, &
                      -0.03963_wp, 0.098224_wp, 0.078945_wp, 0.096585_wp, 0.078472_wp, &
                      0.034254_wp, 0.046298_wp, 0.034521_wp, 0.045395_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.344009_wp, 0.182972_wp, -0.042499_wp, -0.042503_wp, -0.033392_wp, &
                      -0.033381_wp, 0.049425_wp, 0.044244_wp, 0.049479_wp, 0.044175_wp, &
                      0.033125_wp, 0.029611_wp, 0.033137_wp, 0.029615_wp]
         case ("resp")
            charges = [ &
                      -0.639164_wp, 0.617813_wp, -0.20402_wp, -0.200703_wp, -0.038131_wp, &
                      -0.038778_wp, 0.096685_wp, 0.076417_wp, 0.095132_wp, 0.076012_wp, &
                      0.033973_wp, 0.04567_wp, 0.0343_wp, 0.044794_wp]
         case ("mbis")
            charges = [ &
                      -0.610805_wp, 0.677692_wp, -0.431693_wp, -0.431672_wp, -0.188422_wp, &
                      -0.188519_wp, 0.176199_wp, 0.169862_wp, 0.176273_wp, 0.169788_wp, &
                      0.121266_wp, 0.119339_wp, 0.121305_wp, 0.119386_wp]
         case ("chelpg")
            charges = [ &
                      -0.640846_wp, 0.624815_wp, -0.212687_wp, -0.20924_wp, -0.038121_wp, &
                      -0.03864_wp, 0.098951_wp, 0.078768_wp, 0.097349_wp, 0.07832_wp, &
                      0.034561_wp, 0.04642_wp, 0.034842_wp, 0.045508_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (95)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.016932_wp, -0.010044_wp, -0.041457_wp, -0.038836_wp, -0.039473_wp, &
                      -0.039758_wp, -0.010084_wp, 0.016952_wp, -0.041472_wp, -0.038842_wp, &
                      -0.039469_wp, -0.039746_wp, 0.01767_wp, 0.019904_wp, 0.02055_wp, &
                      0.017516_wp, 0.021175_wp, 0.018217_wp, 0.018011_wp, 0.019605_wp, &
                      0.017679_wp, 0.019912_wp, 0.020551_wp, 0.017509_wp, 0.018196_wp, &
                      0.021197_wp, 0.019593_wp, 0.018008_wp]
         case ("resp")
            charges = [ &
                      -0.045173_wp, 0.203821_wp, -0.047148_wp, -0.118581_wp, 0.065839_wp, &
                      -0.259884_wp, 0.220612_wp, -0.049644_wp, -0.060715_wp, -0.099865_wp, &
                      0.070287_wp, -0.270738_wp, 0.036628_wp, 0.005217_wp, 0.021623_wp, &
                      0.038043_wp, 0.009268_wp, -0.014816_wp, 0.049422_wp, 0.051525_wp, &
                      0.039872_wp, 0.00594_wp, 0.015632_wp, 0.033103_wp, -0.015982_wp, &
                      0.006776_wp, 0.054534_wp, 0.054405_wp]
         case ("mbis")
            charges = [ &
                      0.090124_wp, -0.048273_wp, -0.216402_wp, -0.216303_wp, -0.1938_wp, &
                      -0.241036_wp, -0.048178_wp, 0.090141_wp, -0.21638_wp, -0.216538_wp, &
                      -0.193747_wp, -0.240985_wp, 0.105922_wp, 0.100742_wp, 0.104476_wp, &
                      0.101665_wp, 0.102164_wp, 0.101058_wp, 0.103577_wp, 0.105977_wp, &
                      0.10605_wp, 0.100736_wp, 0.104543_wp, 0.101745_wp, 0.101014_wp, &
                      0.102195_wp, 0.105954_wp, 0.103552_wp]
         case ("chelpg")
            charges = [ &
                      -0.04848_wp, 0.217229_wp, -0.046048_wp, -0.133304_wp, 0.07978_wp, &
                      -0.283827_wp, 0.236132_wp, -0.053219_wp, -0.061925_wp, -0.113384_wp, &
                      0.083482_wp, -0.294072_wp, 0.037017_wp, 0.005349_wp, 0.024949_wp, &
                      0.041283_wp, 0.00806_wp, -0.017706_wp, 0.054258_wp, 0.056474_wp, &
                      0.040819_wp, 0.006509_wp, 0.018836_wp, 0.036156_wp, -0.018688_wp, &
                      0.005614_wp, 0.059493_wp, 0.059214_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.017409_wp, -0.00998_wp, -0.041971_wp, -0.039535_wp, -0.040196_wp, &
                      -0.040433_wp, -0.01002_wp, 0.017429_wp, -0.041992_wp, -0.03954_wp, &
                      -0.040195_wp, -0.040424_wp, 0.018738_wp, 0.019957_wp, 0.020501_wp, &
                      0.017743_wp, 0.021077_wp, 0.018771_wp, 0.018294_wp, 0.019637_wp, &
                      0.018766_wp, 0.019958_wp, 0.020497_wp, 0.017735_wp, 0.018769_wp, &
                      0.021095_wp, 0.01962_wp, 0.018287_wp]
         case ("resp")
            charges = [ &
                      -0.046805_wp, 0.214147_wp, -0.049869_wp, -0.122233_wp, 0.070784_wp, &
                      -0.274062_wp, 0.231016_wp, -0.05136_wp, -0.063421_wp, -0.103352_wp, &
                      0.075666_wp, -0.285243_wp, 0.039763_wp, 0.004761_wp, 0.021184_wp, &
                      0.039324_wp, 0.00817_wp, -0.015527_wp, 0.052312_wp, 0.05385_wp, &
                      0.042985_wp, 0.005468_wp, 0.015077_wp, 0.034331_wp, -0.016784_wp, &
                      0.005543_wp, 0.056935_wp, 0.05734_wp]
         case ("mbis")
            charges = [ &
                      0.090925_wp, -0.046904_wp, -0.219834_wp, -0.218167_wp, -0.19604_wp, &
                      -0.244446_wp, -0.046791_wp, 0.09094_wp, -0.219863_wp, -0.218378_wp, &
                      -0.19603_wp, -0.244382_wp, 0.108624_wp, 0.101408_wp, 0.104907_wp, &
                      0.102631_wp, 0.102574_wp, 0.102644_wp, 0.104802_wp, 0.106759_wp, &
                      0.1088_wp, 0.101397_wp, 0.104964_wp, 0.102708_wp, 0.102646_wp, &
                      0.102604_wp, 0.106727_wp, 0.104772_wp]
         case ("chelpg")
            charges = [ &
                      -0.050136_wp, 0.227777_wp, -0.048941_wp, -0.137143_wp, 0.085065_wp, &
                      -0.29835_wp, 0.246662_wp, -0.054947_wp, -0.064702_wp, -0.117141_wp, &
                      0.089275_wp, -0.30897_wp, 0.040215_wp, 0.00495_wp, 0.024554_wp, &
                      0.042608_wp, 0.006898_wp, -0.018496_wp, 0.057213_wp, 0.058872_wp, &
                      0.04397_wp, 0.006078_wp, 0.018339_wp, 0.037442_wp, -0.019586_wp, &
                      0.004304_wp, 0.061971_wp, 0.062219_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.017876_wp, -0.009879_wp, -0.042391_wp, -0.040141_wp, -0.040833_wp, &
                      -0.041019_wp, -0.009919_wp, 0.017894_wp, -0.042418_wp, -0.040147_wp, &
                      -0.040837_wp, -0.041014_wp, 0.019865_wp, 0.019924_wp, 0.02037_wp, &
                      0.017907_wp, 0.020836_wp, 0.019366_wp, 0.018544_wp, 0.019584_wp, &
                      0.019914_wp, 0.01992_wp, 0.020362_wp, 0.017896_wp, 0.019383_wp, &
                      0.020853_wp, 0.019565_wp, 0.018534_wp]
         case ("resp")
            charges = [ &
                      -0.048329_wp, 0.223843_wp, -0.052098_wp, -0.125309_wp, 0.075559_wp, &
                      -0.287014_wp, 0.240749_wp, -0.052957_wp, -0.065604_wp, -0.106322_wp, &
                      0.08087_wp, -0.298523_wp, 0.042844_wp, 0.00406_wp, 0.020549_wp, &
                      0.040302_wp, 0.006795_wp, -0.016081_wp, 0.054925_wp, 0.055787_wp, &
                      0.046039_wp, 0.004748_wp, 0.014344_wp, 0.035265_wp, -0.017429_wp, &
                      0.004042_wp, 0.058946_wp, 0.059998_wp]
         case ("mbis")
            charges = [ &
                      0.0917_wp, -0.045645_wp, -0.222914_wp, -0.219617_wp, -0.198015_wp, &
                      -0.24736_wp, -0.045512_wp, 0.091712_wp, -0.222994_wp, -0.219815_wp, &
                      -0.198045_wp, -0.247286_wp, 0.111359_wp, 0.101825_wp, 0.105094_wp, &
                      0.103395_wp, 0.102644_wp, 0.104247_wp, 0.105878_wp, 0.107285_wp, &
                      0.111584_wp, 0.10181_wp, 0.105144_wp, 0.103467_wp, 0.104291_wp, &
                      0.102676_wp, 0.107246_wp, 0.105841_wp]
         case ("chelpg")
            charges = [ &
                      -0.051679_wp, 0.237657_wp, -0.051306_wp, -0.140377_wp, 0.090143_wp, &
                      -0.311598_wp, 0.256498_wp, -0.056553_wp, -0.066935_wp, -0.120339_wp, &
                      0.094848_wp, -0.322587_wp, 0.043347_wp, 0.004295_wp, 0.023954_wp, &
                      0.04362_wp, 0.005464_wp, -0.019121_wp, 0.059882_wp, 0.060871_wp, &
                      0.047053_wp, 0.00539_wp, 0.017652_wp, 0.038423_wp, -0.020317_wp, &
                      0.002733_wp, 0.064048_wp, 0.064936_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (96)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070371_wp, -0.034197_wp, -0.036422_wp, -0.036255_wp, -0.036205_wp, &
                      -0.036205_wp, -0.036256_wp, -0.036421_wp, -0.034197_wp, -0.070371_wp, &
                      0.021428_wp, 0.021429_wp, 0.023252_wp, 0.019085_wp, 0.019091_wp, &
                      0.018194_wp, 0.018193_wp, 0.018213_wp, 0.018212_wp, 0.018177_wp, &
                      0.018185_wp, 0.018188_wp, 0.018176_wp, 0.018206_wp, 0.01822_wp, &
                      0.018202_wp, 0.018185_wp, 0.019083_wp, 0.019091_wp, 0.023252_wp, &
                      0.021433_wp, 0.021423_wp]
         case ("resp")
            charges = [ &
                      -0.360305_wp, 0.211924_wp, 0.02336_wp, -0.086436_wp, 0.086524_wp, &
                      0.094195_wp, -0.096234_wp, 0.02764_wp, 0.203203_wp, -0.352588_wp, &
                      0.076829_wp, 0.077351_wp, 0.084216_wp, -0.033812_wp, -0.03335_wp, &
                      -0.005246_wp, -0.005174_wp, 0.004087_wp, 0.004071_wp, -0.022197_wp, &
                      -0.022845_wp, -0.024_wp, -0.023378_wp, 0.006699_wp, 0.005732_wp, &
                      -0.005295_wp, -0.005547_wp, -0.03173_wp, -0.030978_wp, 0.08239_wp, &
                      0.074886_wp, 0.076008_wp]
         case ("mbis")
            charges = [ &
                      -0.421709_wp, -0.160121_wp, -0.211213_wp, -0.214054_wp, -0.208478_wp, &
                      -0.208476_wp, -0.21405_wp, -0.211213_wp, -0.160122_wp, -0.421717_wp, &
                      0.127483_wp, 0.127496_wp, 0.129661_wp, 0.099833_wp, 0.099835_wp, &
                      0.105535_wp, 0.105565_wp, 0.105288_wp, 0.105295_wp, 0.104793_wp, &
                      0.104802_wp, 0.104807_wp, 0.104788_wp, 0.105272_wp, 0.105311_wp, &
                      0.10555_wp, 0.105549_wp, 0.099837_wp, 0.099827_wp, 0.129665_wp, &
                      0.127511_wp, 0.12747_wp]
         case ("chelpg")
            charges = [ &
                      -0.37951_wp, 0.222777_wp, 0.030459_wp, -0.10338_wp, 0.09672_wp, &
                      0.105062_wp, -0.113885_wp, 0.034983_wp, 0.213958_wp, -0.371749_wp, &
                      0.08103_wp, 0.081509_wp, 0.088964_wp, -0.03571_wp, -0.035209_wp, &
                      -0.00664_wp, -0.006586_wp, 0.006808_wp, 0.006746_wp, -0.024224_wp, &
                      -0.024882_wp, -0.026148_wp, -0.025491_wp, 0.009568_wp, 0.008621_wp, &
                      -0.006737_wp, -0.007005_wp, -0.033582_wp, -0.032848_wp, 0.087128_wp, &
                      0.07907_wp, 0.080183_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071716_wp, -0.034669_wp, -0.036896_wp, -0.036765_wp, -0.036731_wp, &
                      -0.036731_wp, -0.036765_wp, -0.036896_wp, -0.034669_wp, -0.071716_wp, &
                      0.021695_wp, 0.021696_wp, 0.02349_wp, 0.019439_wp, 0.019444_wp, &
                      0.018536_wp, 0.018534_wp, 0.018518_wp, 0.018517_wp, 0.018455_wp, &
                      0.018461_wp, 0.018463_wp, 0.018453_wp, 0.018513_wp, 0.018523_wp, &
                      0.018541_wp, 0.018528_wp, 0.019437_wp, 0.019444_wp, 0.023491_wp, &
                      0.021699_wp, 0.021692_wp]
         case ("resp")
            charges = [ &
                      -0.371382_wp, 0.221682_wp, 0.02636_wp, -0.086584_wp, 0.092208_wp, &
                      0.100546_wp, -0.096588_wp, 0.03064_wp, 0.213299_wp, -0.363788_wp, &
                      0.078786_wp, 0.079304_wp, 0.086321_wp, -0.035913_wp, -0.035422_wp, &
                      -0.006473_wp, -0.006407_wp, 0.002941_wp, 0.002908_wp, -0.024395_wp, &
                      -0.025076_wp, -0.026327_wp, -0.025677_wp, 0.005586_wp, 0.00457_wp, &
                      -0.006591_wp, -0.006827_wp, -0.033913_wp, -0.03315_wp, 0.084516_wp, &
                      0.076854_wp, 0.077992_wp]
         case ("mbis")
            charges = [ &
                      -0.428955_wp, -0.160408_wp, -0.213207_wp, -0.216072_wp, -0.210383_wp, &
                      -0.210381_wp, -0.216069_wp, -0.213208_wp, -0.160409_wp, -0.428963_wp, &
                      0.129338_wp, 0.129349_wp, 0.131473_wp, 0.100738_wp, 0.100739_wp, &
                      0.106615_wp, 0.106642_wp, 0.106306_wp, 0.10631_wp, 0.10576_wp, &
                      0.105766_wp, 0.10577_wp, 0.105755_wp, 0.106293_wp, 0.106323_wp, &
                      0.106625_wp, 0.106631_wp, 0.100742_wp, 0.100732_wp, 0.131477_wp, &
                      0.129362_wp, 0.129327_wp]
         case ("chelpg")
            charges = [ &
                      -0.390751_wp, 0.232544_wp, 0.033853_wp, -0.103891_wp, 0.102686_wp, &
                      0.111699_wp, -0.114602_wp, 0.038354_wp, 0.224085_wp, -0.383127_wp, &
                      0.083019_wp, 0.083493_wp, 0.091118_wp, -0.037821_wp, -0.037289_wp, &
                      -0.007959_wp, -0.007911_wp, 0.005701_wp, 0.005623_wp, -0.026497_wp, &
                      -0.027187_wp, -0.02855_wp, -0.027863_wp, 0.008495_wp, 0.007499_wp, &
                      -0.00812_wp, -0.008372_wp, -0.035779_wp, -0.035033_wp, 0.089305_wp, &
                      0.081072_wp, 0.082203_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073158_wp, -0.035139_wp, -0.037371_wp, -0.037278_wp, -0.037273_wp, &
                      -0.037273_wp, -0.037279_wp, -0.037371_wp, -0.03514_wp, -0.073158_wp, &
                      0.021952_wp, 0.021951_wp, 0.023687_wp, 0.01983_wp, 0.019834_wp, &
                      0.018902_wp, 0.018898_wp, 0.018847_wp, 0.018844_wp, 0.018741_wp, &
                      0.018744_wp, 0.018746_wp, 0.01874_wp, 0.018844_wp, 0.018847_wp, &
                      0.018904_wp, 0.018896_wp, 0.019828_wp, 0.019834_wp, 0.023687_wp, &
                      0.021952_wp, 0.021951_wp]
         case ("resp")
            charges = [ &
                      -0.383263_wp, 0.232946_wp, 0.030109_wp, -0.086157_wp, 0.098978_wp, &
                      0.108077_wp, -0.0964_wp, 0.034385_wp, 0.224944_wp, -0.375813_wp, &
                      0.080779_wp, 0.081293_wp, 0.08844_wp, -0.03841_wp, -0.037885_wp, &
                      -0.008046_wp, -0.007988_wp, 0.001424_wp, 0.001374_wp, -0.027086_wp, &
                      -0.027805_wp, -0.029168_wp, -0.028482_wp, 0.00411_wp, 0.003034_wp, &
                      -0.008243_wp, -0.008458_wp, -0.0365_wp, -0.035725_wp, 0.086659_wp, &
                      0.078858_wp, 0.08002_wp]
         case ("mbis")
            charges = [ &
                      -0.436646_wp, -0.160624_wp, -0.215257_wp, -0.218152_wp, -0.212355_wp, &
                      -0.212353_wp, -0.218149_wp, -0.215259_wp, -0.160624_wp, -0.436654_wp, &
                      0.131268_wp, 0.131278_wp, 0.133305_wp, 0.101705_wp, 0.101704_wp, &
                      0.107747_wp, 0.10777_wp, 0.107374_wp, 0.107376_wp, 0.106758_wp, &
                      0.106758_wp, 0.106762_wp, 0.106754_wp, 0.107366_wp, 0.107384_wp, &
                      0.10775_wp, 0.107767_wp, 0.10171_wp, 0.101696_wp, 0.133309_wp, &
                      0.131287_wp, 0.131261_wp]
         case ("chelpg")
            charges = [ &
                      -0.402797_wp, 0.243804_wp, 0.03805_wp, -0.103846_wp, 0.10976_wp, &
                      0.119535_wp, -0.114796_wp, 0.042522_wp, 0.235748_wp, -0.395329_wp, &
                      0.085044_wp, 0.085513_wp, 0.093289_wp, -0.040326_wp, -0.039761_wp, &
                      -0.009639_wp, -0.0096_wp, 0.004223_wp, 0.004126_wp, -0.029268_wp, &
                      -0.029998_wp, -0.031471_wp, -0.030747_wp, 0.00706_wp, 0.006004_wp, &
                      -0.009874_wp, -0.010105_wp, -0.038378_wp, -0.037621_wp, 0.091501_wp, &
                      0.083111_wp, 0.084264_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (97)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.035661_wp, -0.074821_wp, -0.074822_wp, 0.056995_wp, 0.056995_wp]
         case ("resp")
            charges = [ &
                      -0.372967_wp, -0.046845_wp, -0.045048_wp, 0.232341_wp, 0.232519_wp]
         case ("mbis")
            charges = [ &
                      -0.228872_wp, -0.090617_wp, -0.090619_wp, 0.205058_wp, 0.205058_wp]
         case ("chelpg")
            charges = [ &
                      -0.377288_wp, -0.046091_wp, -0.044299_wp, 0.233752_wp, 0.233925_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.045471_wp, -0.092911_wp, -0.092912_wp, 0.07018_wp, 0.07018_wp]
         case ("resp")
            charges = [ &
                      -0.36195_wp, -0.079835_wp, -0.077927_wp, 0.259787_wp, 0.259925_wp]
         case ("mbis")
            charges = [ &
                      -0.212269_wp, -0.117295_wp, -0.117298_wp, 0.223435_wp, 0.223435_wp]
         case ("chelpg")
            charges = [ &
                      -0.36607_wp, -0.079121_wp, -0.077217_wp, 0.261137_wp, 0.26127_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.047661_wp, -0.096676_wp, -0.096678_wp, 0.07285_wp, 0.07285_wp]
         case ("resp")
            charges = [ &
                      -0.358573_wp, -0.087188_wp, -0.085252_wp, 0.265442_wp, 0.26557_wp]
         case ("mbis")
            charges = [ &
                      -0.208205_wp, -0.122896_wp, -0.122899_wp, 0.227004_wp, 0.227004_wp]
         case ("chelpg")
            charges = [ &
                      -0.362653_wp, -0.086483_wp, -0.084551_wp, 0.266782_wp, 0.266905_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (98)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.190625_wp, 0.04103_wp, -0.038268_wp, -0.032475_wp, -0.068644_wp, &
                      0.04103_wp, -0.038268_wp, -0.032475_wp, -0.068644_wp, 0.023123_wp, &
                      0.024866_wp, 0.024866_wp, 0.024866_wp, 0.024865_wp, 0.016898_wp, &
                      0.016898_wp, 0.016895_wp, 0.0169_wp, 0.023123_wp, 0.023122_wp, &
                      0.023954_wp, 0.019975_wp, 0.019975_wp, 0.019974_wp, 0.019976_wp, &
                      0.023122_wp, 0.023954_wp]
         case ("resp")
            charges = [ &
                      -0.436058_wp, 0.16782_wp, -0.03873_wp, 0.176665_wp, -0.275681_wp, &
                      0.165625_wp, -0.038593_wp, 0.176018_wp, -0.280117_wp, 0.061005_wp, &
                      0.030222_wp, 0.030306_wp, 0.030751_wp, 0.030649_wp, 0.004735_wp, &
                      0.004003_wp, 0.004201_wp, 0.005315_wp, 0.059728_wp, 0.059815_wp, &
                      0.066927_wp, -0.033665_wp, -0.033827_wp, -0.033079_wp, -0.033147_wp, &
                      0.060871_wp, 0.068242_wp]
         case ("mbis")
            charges = [ &
                      -0.364424_wp, 0.079553_wp, -0.286151_wp, -0.148539_wp, -0.420827_wp, &
                      0.079551_wp, -0.28615_wp, -0.148538_wp, -0.420828_wp, 0.129813_wp, &
                      0.128562_wp, 0.128562_wp, 0.128556_wp, 0.128566_wp, 0.056644_wp, &
                      0.056645_wp, 0.056638_wp, 0.056648_wp, 0.129813_wp, 0.129806_wp, &
                      0.129902_wp, 0.099127_wp, 0.099124_wp, 0.099124_wp, 0.099127_wp, &
                      0.129807_wp, 0.129902_wp]
         case ("chelpg")
            charges = [ &
                      -0.438688_wp, 0.176771_wp, -0.049919_wp, 0.190486_wp, -0.29042_wp, &
                      0.174496_wp, -0.049726_wp, 0.1898_wp, -0.294919_wp, 0.064425_wp, &
                      0.032217_wp, 0.032299_wp, 0.032736_wp, 0.032643_wp, 0.002493_wp, &
                      0.001757_wp, 0.001968_wp, 0.003093_wp, 0.063154_wp, 0.063214_wp, &
                      0.070078_wp, -0.036153_wp, -0.036336_wp, -0.03557_wp, -0.035621_wp, &
                      0.064311_wp, 0.071411_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.204652_wp, 0.039807_wp, -0.040294_wp, -0.031753_wp, -0.070175_wp, &
                      0.039807_wp, -0.040295_wp, -0.031753_wp, -0.070175_wp, 0.023513_wp, &
                      0.023564_wp, 0.023565_wp, 0.023564_wp, 0.023562_wp, 0.020947_wp, &
                      0.020949_wp, 0.020946_wp, 0.020949_wp, 0.023513_wp, 0.023512_wp, &
                      0.024695_wp, 0.022002_wp, 0.022001_wp, 0.022001_wp, 0.022002_wp, &
                      0.023512_wp, 0.024695_wp]
         case ("resp")
            charges = [ &
                      -0.483704_wp, 0.179173_wp, -0.045737_wp, 0.201919_wp, -0.292257_wp, &
                      0.176873_wp, -0.045382_wp, 0.201226_wp, -0.296759_wp, 0.064122_wp, &
                      0.028051_wp, 0.028157_wp, 0.028565_wp, 0.028411_wp, 0.011078_wp, &
                      0.010386_wp, 0.010612_wp, 0.011646_wp, 0.06285_wp, 0.062924_wp, &
                      0.069636_wp, -0.036887_wp, -0.037142_wp, -0.03639_wp, -0.036365_wp, &
                      0.064016_wp, 0.070977_wp]
         case ("mbis")
            charges = [ &
                      -0.387999_wp, 0.073272_wp, -0.291945_wp, -0.145085_wp, -0.43121_wp, &
                      0.07327_wp, -0.291944_wp, -0.145084_wp, -0.431211_wp, 0.132293_wp, &
                      0.128013_wp, 0.128015_wp, 0.128009_wp, 0.128016_wp, 0.065594_wp, &
                      0.065597_wp, 0.065591_wp, 0.065597_wp, 0.132293_wp, 0.132288_wp, &
                      0.132934_wp, 0.102123_wp, 0.102119_wp, 0.10212_wp, 0.102121_wp, &
                      0.132288_wp, 0.132934_wp]
         case ("chelpg")
            charges = [ &
                      -0.486186_wp, 0.188073_wp, -0.057402_wp, 0.216315_wp, -0.307374_wp, &
                      0.185682_wp, -0.056971_wp, 0.215576_wp, -0.311937_wp, 0.067621_wp, &
                      0.030145_wp, 0.030248_wp, 0.030642_wp, 0.030498_wp, 0.008864_wp, &
                      0.008169_wp, 0.008409_wp, 0.009453_wp, 0.066357_wp, 0.066402_wp, &
                      0.072851_wp, -0.039484_wp, -0.039761_wp, -0.03899_wp, -0.038946_wp, &
                      0.067536_wp, 0.074209_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.21324_wp, 0.039187_wp, -0.041593_wp, -0.031138_wp, -0.07097_wp, &
                      0.039187_wp, -0.041594_wp, -0.031139_wp, -0.07097_wp, 0.023541_wp, &
                      0.022277_wp, 0.02228_wp, 0.022279_wp, 0.022276_wp, 0.023648_wp, &
                      0.02365_wp, 0.023648_wp, 0.023649_wp, 0.023541_wp, 0.023541_wp, &
                      0.025243_wp, 0.023482_wp, 0.023481_wp, 0.023481_wp, 0.023481_wp, &
                      0.023541_wp, 0.025243_wp]
         case ("resp")
            charges = [ &
                      -0.513849_wp, 0.187653_wp, -0.051762_wp, 0.217807_wp, -0.301772_wp, &
                      0.18528_wp, -0.051271_wp, 0.217091_wp, -0.306314_wp, 0.065592_wp, &
                      0.026237_wp, 0.02636_wp, 0.026744_wp, 0.026554_wp, 0.015165_wp, &
                      0.014504_wp, 0.014751_wp, 0.015727_wp, 0.064326_wp, 0.064388_wp, &
                      0.071407_wp, -0.038384_wp, -0.038696_wp, -0.037944_wp, -0.037861_wp, &
                      0.065504_wp, 0.072763_wp]
         case ("mbis")
            charges = [ &
                      -0.4028_wp, 0.069694_wp, -0.294908_wp, -0.143075_wp, -0.43655_wp, &
                      0.069691_wp, -0.294906_wp, -0.143074_wp, -0.436551_wp, 0.133262_wp, &
                      0.126757_wp, 0.126762_wp, 0.126755_wp, 0.126761_wp, 0.071372_wp, &
                      0.071377_wp, 0.071372_wp, 0.071373_wp, 0.133262_wp, 0.133258_wp, &
                      0.134792_wp, 0.104337_wp, 0.104332_wp, 0.104334_wp, 0.104334_wp, &
                      0.133258_wp, 0.134792_wp]
         case ("chelpg")
            charges = [ &
                      -0.516237_wp, 0.196548_wp, -0.063761_wp, 0.232535_wp, -0.317088_wp, &
                      0.194079_wp, -0.063182_wp, 0.231768_wp, -0.321691_wp, 0.069133_wp, &
                      0.028403_wp, 0.028523_wp, 0.028891_wp, 0.028711_wp, 0.012962_wp, &
                      0.012298_wp, 0.012559_wp, 0.013546_wp, 0.067875_wp, 0.067909_wp, &
                      0.074651_wp, -0.041041_wp, -0.041376_wp, -0.040604_wp, -0.040502_wp, &
                      0.069067_wp, 0.076025_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (99)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.033842_wp, 0.033842_wp, -0.03049_wp, -0.030491_wp, -0.027765_wp, &
                      -0.027765_wp, -0.067125_wp, -0.067124_wp, 0.047514_wp, 0.047514_wp, &
                      0.044024_wp, 0.044023_wp]
         case ("resp")
            charges = [ &
                      0.075072_wp, 0.068903_wp, -0.115249_wp, -0.116009_wp, -0.140377_wp, &
                      -0.13994_wp, -0.092486_wp, -0.090648_wp, 0.131828_wp, 0.133135_wp, &
                      0.142815_wp, 0.142956_wp]
         case ("mbis")
            charges = [ &
                      0.101628_wp, 0.101629_wp, -0.159745_wp, -0.159746_wp, -0.133206_wp, &
                      -0.133203_wp, -0.123049_wp, -0.123049_wp, 0.167975_wp, 0.167975_wp, &
                      0.146396_wp, 0.146395_wp]
         case ("chelpg")
            charges = [ &
                      0.076958_wp, 0.070389_wp, -0.118165_wp, -0.118537_wp, -0.139619_wp, &
                      -0.139519_wp, -0.093024_wp, -0.091131_wp, 0.132779_wp, 0.134012_wp, &
                      0.142829_wp, 0.143027_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.027875_wp, 0.027875_wp, -0.030889_wp, -0.03089_wp, -0.025265_wp, &
                      -0.025265_wp, -0.082617_wp, -0.082616_wp, 0.056096_wp, 0.056096_wp, &
                      0.054798_wp, 0.054798_wp]
         case ("resp")
            charges = [ &
                      0.070529_wp, 0.064496_wp, -0.12062_wp, -0.121391_wp, -0.150806_wp, &
                      -0.150185_wp, -0.114234_wp, -0.112459_wp, 0.152168_wp, 0.153411_wp, &
                      0.1645_wp, 0.164591_wp]
         case ("mbis")
            charges = [ &
                      0.096029_wp, 0.09603_wp, -0.168637_wp, -0.168637_wp, -0.138058_wp, &
                      -0.138056_wp, -0.140134_wp, -0.140134_wp, 0.184949_wp, 0.18495_wp, &
                      0.16585_wp, 0.165848_wp]
         case ("chelpg")
            charges = [ &
                      0.072444_wp, 0.065997_wp, -0.123569_wp, -0.123937_wp, -0.150045_wp, &
                      -0.149774_wp, -0.114782_wp, -0.112949_wp, 0.153132_wp, 0.154297_wp, &
                      0.164518_wp, 0.164668_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.026905_wp, 0.026904_wp, -0.030874_wp, -0.030875_wp, -0.024605_wp, &
                      -0.024605_wp, -0.085011_wp, -0.08501_wp, 0.057215_wp, 0.057215_wp, &
                      0.05637_wp, 0.05637_wp]
         case ("resp")
            charges = [ &
                      0.069647_wp, 0.063641_wp, -0.121184_wp, -0.12196_wp, -0.151678_wp, &
                      -0.151023_wp, -0.117603_wp, -0.115839_wp, 0.154898_wp, 0.15613_wp, &
                      0.167445_wp, 0.167527_wp]
         case ("mbis")
            charges = [ &
                      0.094966_wp, 0.094967_wp, -0.169691_wp, -0.169692_wp, -0.13814_wp, &
                      -0.138138_wp, -0.142735_wp, -0.142735_wp, 0.187113_wp, 0.187114_wp, &
                      0.168485_wp, 0.168484_wp]
         case ("chelpg")
            charges = [ &
                      0.071567_wp, 0.065144_wp, -0.124137_wp, -0.124508_wp, -0.150913_wp, &
                      -0.150609_wp, -0.118156_wp, -0.116335_wp, 0.155862_wp, 0.157016_wp, &
                      0.167464_wp, 0.167605_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (100)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.004729_wp, 0.004729_wp, 0.061659_wp, -0.066389_wp, -0.066388_wp, &
                      0.06166_wp]
         case ("resp")
            charges = [ &
                      -0.125968_wp, -0.132591_wp, 0.199396_wp, -0.071786_wp, -0.07019_wp, &
                      0.201138_wp]
         case ("mbis")
            charges = [ &
                      -0.083087_wp, -0.083088_wp, 0.194011_wp, -0.110924_wp, -0.110922_wp, &
                      0.194011_wp]
         case ("chelpg")
            charges = [ &
                      -0.126335_wp, -0.133009_wp, 0.199689_wp, -0.071698_wp, -0.070091_wp, &
                      0.201445_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.005291_wp, 0.005291_wp, 0.067554_wp, -0.072845_wp, -0.072844_wp, &
                      0.067554_wp]
         case ("resp")
            charges = [ &
                      -0.127631_wp, -0.134343_wp, 0.211283_wp, -0.081989_wp, -0.080376_wp, &
                      0.213057_wp]
         case ("mbis")
            charges = [ &
                      -0.085972_wp, -0.085973_wp, 0.204819_wp, -0.118847_wp, -0.118846_wp, &
                      0.204819_wp]
         case ("chelpg")
            charges = [ &
                      -0.127996_wp, -0.134757_wp, 0.211576_wp, -0.081905_wp, -0.080281_wp, &
                      0.213363_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.006142_wp, 0.006142_wp, 0.073378_wp, -0.07952_wp, -0.079519_wp, &
                      0.073378_wp]
         case ("resp")
            charges = [ &
                      -0.127962_wp, -0.13476_wp, 0.222736_wp, -0.093092_wp, -0.091464_wp, &
                      0.224542_wp]
         case ("mbis")
            charges = [ &
                      -0.088087_wp, -0.088088_wp, 0.215279_wp, -0.127192_wp, -0.12719_wp, &
                      0.215279_wp]
         case ("chelpg")
            charges = [ &
                      -0.128323_wp, -0.13517_wp, 0.223029_wp, -0.093011_wp, -0.091372_wp, &
                      0.224847_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (101)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      1.2e-05_wp, 1.1e-05_wp, -0.058913_wp, 0.058901_wp, -0.058912_wp, &
                      0.058901_wp]
         case ("resp")
            charges = [ &
                      -0.125703_wp, -0.13763_wp, -0.043067_wp, 0.171367_wp, -0.038863_wp, &
                      0.173896_wp]
         case ("mbis")
            charges = [ &
                      -0.103844_wp, -0.103842_wp, -0.078948_wp, 0.182791_wp, -0.078948_wp, &
                      0.182791_wp]
         case ("chelpg")
            charges = [ &
                      -0.126013_wp, -0.138101_wp, -0.042933_wp, 0.171584_wp, -0.038696_wp, &
                      0.17416_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.003091_wp, 0.003091_wp, -0.076927_wp, 0.073835_wp, -0.076926_wp, &
                      0.073835_wp]
         case ("resp")
            charges = [ &
                      -0.12937_wp, -0.141286_wp, -0.06893_wp, 0.200929_wp, -0.064817_wp, &
                      0.203473_wp]
         case ("mbis")
            charges = [ &
                      -0.106885_wp, -0.106882_wp, -0.101632_wp, 0.208516_wp, -0.101632_wp, &
                      0.208516_wp]
         case ("chelpg")
            charges = [ &
                      -0.129672_wp, -0.141732_wp, -0.068803_wp, 0.201141_wp, -0.06466_wp, &
                      0.203727_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.003724_wp, 0.003724_wp, -0.079837_wp, 0.076113_wp, -0.079836_wp, &
                      0.076113_wp]
         case ("resp")
            charges = [ &
                      -0.129409_wp, -0.141318_wp, -0.073326_wp, 0.20537_wp, -0.069233_wp, &
                      0.207916_wp]
         case ("mbis")
            charges = [ &
                      -0.10692_wp, -0.106917_wp, -0.105381_wp, 0.2123_wp, -0.105381_wp, &
                      0.2123_wp]
         case ("chelpg")
            charges = [ &
                      -0.129709_wp, -0.141761_wp, -0.073201_wp, 0.205581_wp, -0.069078_wp, &
                      0.208168_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (102)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.066792_wp, 0.058084_wp, 0.058081_wp, -0.091476_wp, -0.091481_wp]
         case ("resp")
            charges = [ &
                      -0.301285_wp, 0.218036_wp, 0.217766_wp, -0.067626_wp, -0.066891_wp]
         case ("mbis")
            charges = [ &
                      -0.081185_wp, 0.182812_wp, 0.182811_wp, -0.142218_wp, -0.142221_wp]
         case ("chelpg")
            charges = [ &
                      -0.305708_wp, 0.219441_wp, 0.219166_wp, -0.066815_wp, -0.066084_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.075578_wp, 0.071957_wp, 0.071953_wp, -0.109742_wp, -0.109747_wp]
         case ("resp")
            charges = [ &
                      -0.283383_wp, 0.244017_wp, 0.243738_wp, -0.102609_wp, -0.101763_wp]
         case ("mbis")
            charges = [ &
                      -0.069602_wp, 0.203333_wp, 0.203331_wp, -0.168529_wp, -0.168533_wp]
         case ("chelpg")
            charges = [ &
                      -0.287626_wp, 0.245368_wp, 0.245085_wp, -0.101835_wp, -0.100993_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.077063_wp, 0.07411_wp, 0.074106_wp, -0.112637_wp, -0.112642_wp]
         case ("resp")
            charges = [ &
                      -0.279957_wp, 0.248088_wp, 0.247808_wp, -0.108402_wp, -0.107536_wp]
         case ("mbis")
            charges = [ &
                      -0.067379_wp, 0.206427_wp, 0.206424_wp, -0.172734_wp, -0.172738_wp]
         case ("chelpg")
            charges = [ &
                      -0.284175_wp, 0.249433_wp, 0.249148_wp, -0.107635_wp, -0.106771_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (103)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070115_wp, 0.043238_wp, -0.192419_wp, 0.043244_wp, -0.070115_wp, &
                      0.028394_wp, 0.029206_wp, 0.029199_wp, 0.018143_wp, 0.01814_wp, &
                      0.018143_wp, 0.018144_wp, 0.028394_wp, 0.029202_wp, 0.029203_wp]
         case ("resp")
            charges = [ &
                      -0.335436_wp, 0.308399_wp, -0.422582_wp, 0.307559_wp, -0.333307_wp, &
                      0.076376_wp, 0.098697_wp, 0.098349_wp, -0.017321_wp, -0.017527_wp, &
                      -0.017652_wp, -0.017525_wp, 0.076152_wp, 0.098106_wp, 0.097713_wp]
         case ("mbis")
            charges = [ &
                      -0.489632_wp, 0.124408_wp, -0.36736_wp, 0.124391_wp, -0.489627_wp, &
                      0.142495_wp, 0.150897_wp, 0.150879_wp, 0.052319_wp, 0.052315_wp, &
                      0.052319_wp, 0.052325_wp, 0.142497_wp, 0.150883_wp, 0.150893_wp]
         case ("chelpg")
            charges = [ &
                      -0.349195_wp, 0.321288_wp, -0.427785_wp, 0.320492_wp, -0.347098_wp, &
                      0.079296_wp, 0.1019_wp, 0.101544_wp, -0.020238_wp, -0.020459_wp, &
                      -0.020589_wp, -0.020464_wp, 0.079078_wp, 0.101313_wp, 0.100918_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073036_wp, 0.042271_wp, -0.208753_wp, 0.042277_wp, -0.073036_wp, &
                      0.033381_wp, 0.027792_wp, 0.027791_wp, 0.02309_wp, 0.023082_wp, &
                      0.02309_wp, 0.023087_wp, 0.033381_wp, 0.027793_wp, 0.02779_wp]
         case ("resp")
            charges = [ &
                      -0.354993_wp, 0.323917_wp, -0.477028_wp, 0.322837_wp, -0.352995_wp, &
                      0.089638_wp, 0.099754_wp, 0.099409_wp, -0.009363_wp, -0.009578_wp, &
                      -0.009607_wp, -0.009502_wp, 0.089453_wp, 0.099229_wp, 0.098828_wp]
         case ("mbis")
            charges = [ &
                      -0.504615_wp, 0.120699_wp, -0.395509_wp, 0.12068_wp, -0.504609_wp, &
                      0.154226_wp, 0.151368_wp, 0.15136_wp, 0.062365_wp, 0.062352_wp, &
                      0.062365_wp, 0.062364_wp, 0.154228_wp, 0.151362_wp, 0.151366_wp]
         case ("chelpg")
            charges = [ &
                      -0.36884_wp, 0.336601_wp, -0.482085_wp, 0.335563_wp, -0.366876_wp, &
                      0.092598_wp, 0.102986_wp, 0.102632_wp, -0.012219_wp, -0.012449_wp, &
                      -0.012482_wp, -0.012379_wp, 0.092421_wp, 0.102466_wp, 0.102063_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074092_wp, 0.041996_wp, -0.215068_wp, 0.042002_wp, -0.074092_wp, &
                      0.035542_wp, 0.026908_wp, 0.02691_wp, 0.025138_wp, 0.025127_wp, &
                      0.025138_wp, 0.025133_wp, 0.035543_wp, 0.02691_wp, 0.026907_wp]
         case ("resp")
            charges = [ &
                      -0.362194_wp, 0.330524_wp, -0.499084_wp, 0.329346_wp, -0.360244_wp, &
                      0.095246_wp, 0.099512_wp, 0.099167_wp, -0.006106_wp, -0.006327_wp, &
                      -0.006317_wp, -0.00622_wp, 0.095075_wp, 0.099012_wp, 0.098608_wp]
         case ("mbis")
            charges = [ &
                      -0.509576_wp, 0.119352_wp, -0.406652_wp, 0.119333_wp, -0.50957_wp, &
                      0.159063_wp, 0.150823_wp, 0.150819_wp, 0.066431_wp, 0.066414_wp, &
                      0.066431_wp, 0.066427_wp, 0.159065_wp, 0.15082_wp, 0.150821_wp]
         case ("chelpg")
            charges = [ &
                      -0.376063_wp, 0.343109_wp, -0.504079_wp, 0.341972_wp, -0.374147_wp, &
                      0.098221_wp, 0.102751_wp, 0.102398_wp, -0.008932_wp, -0.009168_wp, &
                      -0.009161_wp, -0.009067_wp, 0.098059_wp, 0.102257_wp, 0.101851_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (104)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.088438_wp, -0.026509_wp, 0.035905_wp, 0.027552_wp, -0.071202_wp, &
                      0.024989_wp, 0.027202_wp, 0.026317_wp, -0.026507_wp, 0.035861_wp, &
                      0.027572_wp, -0.071216_wp, 0.027185_wp, 0.026314_wp, 0.024982_wp]
         case ("resp")
            charges = [ &
                      -0.362149_wp, 0.168135_wp, 0.069759_wp, -0.00242_wp, -0.304458_wp, &
                      0.104764_wp, 0.079708_wp, 0.066522_wp, 0.156138_wp, 0.073015_wp, &
                      0.001875_wp, -0.305489_wp, 0.081074_wp, 0.067437_wp, 0.106091_wp]
         case ("mbis")
            charges = [ &
                      -0.202922_wp, -0.193492_wp, 0.153755_wp, 0.134219_wp, -0.401724_wp, &
                      0.133761_wp, 0.146375_wp, 0.12859_wp, -0.19345_wp, 0.153734_wp, &
                      0.134202_wp, -0.401743_wp, 0.146295_wp, 0.128587_wp, 0.133817_wp]
         case ("chelpg")
            charges = [ &
                      -0.367261_wp, 0.18382_wp, 0.067085_wp, -0.006598_wp, -0.324654_wp, &
                      0.109886_wp, 0.084055_wp, 0.071004_wp, 0.171577_wp, 0.070414_wp, &
                      -0.002226_wp, -0.325704_wp, 0.085469_wp, 0.071937_wp, 0.111197_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.123878_wp, -0.024573_wp, 0.039563_wp, 0.034081_wp, -0.071776_wp, &
                      0.028436_wp, 0.025524_wp, 0.030721_wp, -0.024556_wp, 0.039557_wp, &
                      0.034029_wp, -0.071783_wp, 0.025497_wp, 0.030705_wp, 0.028459_wp]
         case ("resp")
            charges = [ &
                      -0.442984_wp, 0.197479_wp, 0.072163_wp, 0.004229_wp, -0.314714_wp, &
                      0.112281_wp, 0.076995_wp, 0.073977_wp, 0.18552_wp, 0.075463_wp, &
                      0.008451_wp, -0.315697_wp, 0.078333_wp, 0.074848_wp, 0.113656_wp]
         case ("mbis")
            charges = [ &
                      -0.265703_wp, -0.182448_wp, 0.159781_wp, 0.145038_wp, -0.414115_wp, &
                      0.141626_wp, 0.145333_wp, 0.137656_wp, -0.182302_wp, 0.159814_wp, &
                      0.144871_wp, -0.414144_wp, 0.145233_wp, 0.137632_wp, 0.141733_wp]
         case ("chelpg")
            charges = [ &
                      -0.448165_wp, 0.213359_wp, 0.069456_wp, 6e-06_wp, -0.33519_wp, &
                      0.117473_wp, 0.081405_wp, 0.078522_wp, 0.201187_wp, 0.07282_wp, &
                      0.004295_wp, -0.336206_wp, 0.082792_wp, 0.079413_wp, 0.118834_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.133958_wp, -0.02394_wp, 0.040283_wp, 0.036035_wp, -0.071763_wp, &
                      0.029503_wp, 0.02481_wp, 0.032087_wp, -0.023919_wp, 0.040288_wp, &
                      0.035963_wp, -0.071767_wp, 0.024781_wp, 0.032068_wp, 0.029535_wp]
         case ("resp")
            charges = [ &
                      -0.466278_wp, 0.206245_wp, 0.072117_wp, 0.006349_wp, -0.317157_wp, &
                      0.114558_wp, 0.075675_wp, 0.076265_wp, 0.194315_wp, 0.075425_wp, &
                      0.010536_wp, -0.318095_wp, 0.076996_wp, 0.077112_wp, 0.115936_wp]
         case ("mbis")
            charges = [ &
                      -0.28369_wp, -0.179002_wp, 0.160875_wp, 0.148219_wp, -0.416788_wp, &
                      0.143851_wp, 0.144434_wp, 0.140274_wp, -0.178826_wp, 0.160922_wp, &
                      0.14801_wp, -0.41682_wp, 0.144329_wp, 0.140242_wp, 0.143974_wp]
         case ("chelpg")
            charges = [ &
                      -0.47147_wp, 0.222159_wp, 0.069405_wp, 0.002118_wp, -0.33769_wp, &
                      0.119763_wp, 0.080098_wp, 0.080824_wp, 0.210024_wp, 0.072774_wp, &
                      0.006371_wp, -0.338665_wp, 0.081469_wp, 0.081692_wp, 0.121128_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (105)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.183203_wp, 0.010196_wp, 0.010195_wp, 0.00817_wp, 0.024322_wp, &
                      0.088746_wp, 0.024323_wp, 0.008169_wp, -0.070462_wp, 0.026268_wp, &
                      0.022339_wp, 0.026401_wp, -0.070464_wp, 0.02234_wp, 0.026268_wp, &
                      0.026397_wp]
         case ("resp")
            charges = [ &
                      -0.742255_wp, 0.354721_wp, 0.362177_wp, -0.052265_wp, 0.002883_wp, &
                      0.322608_wp, 0.001318_wp, -0.053159_wp, -0.360488_wp, 0.095621_wp, &
                      0.091031_wp, 0.077214_wp, -0.367739_wp, 0.092053_wp, 0.097674_wp, &
                      0.078607_wp]
         case ("mbis")
            charges = [ &
                      -0.599539_wp, 0.076312_wp, 0.07632_wp, 0.040151_wp, 0.082021_wp, &
                      0.317989_wp, 0.082024_wp, 0.040142_wp, -0.479858_wp, 0.14964_wp, &
                      0.131774_wp, 0.140733_wp, -0.479853_wp, 0.131782_wp, 0.149642_wp, &
                      0.140725_wp]
         case ("chelpg")
            charges = [ &
                      -0.748995_wp, 0.368632_wp, 0.376289_wp, -0.055473_wp, -0.00014_wp, &
                      0.323103_wp, -0.001739_wp, -0.056376_wp, -0.374449_wp, 0.098984_wp, &
                      0.094088_wp, 0.08017_wp, -0.381981_wp, 0.09517_wp, 0.1011_wp, &
                      0.081618_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.195701_wp, 0.009133_wp, 0.009133_wp, 0.011626_wp, 0.025034_wp, &
                      0.095648_wp, 0.025036_wp, 0.011626_wp, -0.072727_wp, 0.023047_wp, &
                      0.0253_wp, 0.028619_wp, -0.07273_wp, 0.0253_wp, 0.023047_wp, &
                      0.028615_wp]
         case ("resp")
            charges = [ &
                      -0.837286_wp, 0.39074_wp, 0.398677_wp, -0.049045_wp, 3.4e-05_wp, &
                      0.361062_wp, -0.001602_wp, -0.050051_wp, -0.384987_wp, 0.094509_wp, &
                      0.101383_wp, 0.08431_wp, -0.392636_wp, 0.102485_wp, 0.096629_wp, &
                      0.085778_wp]
         case ("mbis")
            charges = [ &
                      -0.650994_wp, 0.085815_wp, 0.085826_wp, 0.045861_wp, 0.083193_wp, &
                      0.342141_wp, 0.083197_wp, 0.045852_wp, -0.49406_wp, 0.146771_wp, &
                      0.139823_wp, 0.147023_wp, -0.494061_wp, 0.139829_wp, 0.146774_wp, &
                      0.147016_wp]
         case ("chelpg")
            charges = [ &
                      -0.844049_wp, 0.404752_wp, 0.412892_wp, -0.052277_wp, -0.003009_wp, &
                      0.361553_wp, -0.00468_wp, -0.053293_wp, -0.399112_wp, 0.097913_wp, &
                      0.104479_wp, 0.087307_wp, -0.40704_wp, 0.10564_wp, 0.100095_wp, &
                      0.088828_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.201746_wp, 0.008714_wp, 0.008714_wp, 0.013865_wp, 0.025012_wp, &
                      0.098796_wp, 0.025014_wp, 0.013866_wp, -0.07371_wp, 0.020928_wp, &
                      0.026876_wp, 0.029796_wp, -0.073714_wp, 0.026875_wp, 0.020927_wp, &
                      0.029793_wp]
         case ("resp")
            charges = [ &
                      -0.884219_wp, 0.408234_wp, 0.416429_wp, -0.046201_wp, -0.002085_wp, &
                      0.379439_wp, -0.003762_wp, -0.047272_wp, -0.39604_wp, 0.092737_wp, &
                      0.106697_wp, 0.087848_wp, -0.403871_wp, 0.107833_wp, 0.094886_wp, &
                      0.089347_wp]
         case ("mbis")
            charges = [ &
                      -0.676027_wp, 0.090455_wp, 0.090468_wp, 0.049668_wp, 0.083034_wp, &
                      0.353364_wp, 0.083039_wp, 0.049659_wp, -0.499956_wp, 0.144203_wp, &
                      0.143832_wp, 0.150094_wp, -0.499959_wp, 0.143837_wp, 0.144205_wp, &
                      0.150089_wp]
         case ("chelpg")
            charges = [ &
                      -0.891012_wp, 0.422293_wp, 0.430692_wp, -0.049436_wp, -0.00514_wp, &
                      0.379938_wp, -0.006852_wp, -0.050518_wp, -0.410231_wp, 0.096153_wp, &
                      0.10981_wp, 0.090861_wp, -0.418342_wp, 0.111005_wp, 0.098365_wp, &
                      0.092414_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (106)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.053734_wp, 0.001287_wp, -0.053841_wp, 0.05315_wp, 0.053138_wp]
         case ("resp")
            charges = [ &
                      -0.019165_wp, -0.425112_wp, -0.0161_wp, 0.23081_wp, 0.229567_wp]
         case ("mbis")
            charges = [ &
                      -0.014962_wp, -0.413362_wp, -0.015055_wp, 0.221693_wp, 0.221686_wp]
         case ("chelpg")
            charges = [ &
                      -0.018446_wp, -0.429485_wp, -0.015382_wp, 0.232283_wp, 0.23103_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068499_wp, 0.009402_wp, -0.068599_wp, 0.063878_wp, 0.063818_wp]
         case ("resp")
            charges = [ &
                      -0.041185_wp, -0.438315_wp, -0.038176_wp, 0.259539_wp, 0.258137_wp]
         case ("mbis")
            charges = [ &
                      -0.03635_wp, -0.400413_wp, -0.036432_wp, 0.23664_wp, 0.236554_wp]
         case ("chelpg")
            charges = [ &
                      -0.040495_wp, -0.442534_wp, -0.037487_wp, 0.260964_wp, 0.259552_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073028_wp, 0.012106_wp, -0.073125_wp, 0.067062_wp, 0.066986_wp]
         case ("resp")
            charges = [ &
                      -0.04881_wp, -0.440249_wp, -0.045843_wp, 0.26818_wp, 0.266722_wp]
         case ("mbis")
            charges = [ &
                      -0.042992_wp, -0.395524_wp, -0.04307_wp, 0.240849_wp, 0.240736_wp]
         case ("chelpg")
            charges = [ &
                      -0.048129_wp, -0.44442_wp, -0.045163_wp, 0.269589_wp, 0.268123_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (107)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073808_wp, 0.072897_wp, 0.025643_wp, 0.027248_wp, 0.028048_wp, &
                      0.018327_wp, -0.077426_wp, 0.02364_wp, 0.026702_wp, 0.025995_wp, &
                      -0.19456_wp, 0.072914_wp, -0.07381_wp, -0.077423_wp, 0.018346_wp, &
                      0.028048_wp, 0.027246_wp, 0.025631_wp, 0.023642_wp, 0.026706_wp, &
                      0.025999_wp]
         case ("resp")
            charges = [ &
                      -0.551905_wp, 0.721058_wp, 0.119269_wp, 0.126423_wp, 0.149215_wp, &
                      -0.074754_wp, -0.642072_wp, 0.145747_wp, 0.161477_wp, 0.138452_wp, &
                      -0.585942_wp, 0.721454_wp, -0.535719_wp, -0.636485_wp, -0.077708_wp, &
                      0.144703_wp, 0.122556_wp, 0.113779_wp, 0.143121_wp, 0.159909_wp, &
                      0.137422_wp]
         case ("mbis")
            charges = [ &
                      -0.527352_wp, 0.305076_wp, 0.146853_wp, 0.153045_wp, 0.156669_wp, &
                      0.055589_wp, -0.54348_wp, 0.143032_wp, 0.157936_wp, 0.150101_wp, &
                      -0.394962_wp, 0.305053_wp, -0.527331_wp, -0.543473_wp, 0.055628_wp, &
                      0.156658_wp, 0.15304_wp, 0.146826_wp, 0.143037_wp, 0.157947_wp, &
                      0.150109_wp]
         case ("chelpg")
            charges = [ &
                      -0.570649_wp, 0.744258_wp, 0.123185_wp, 0.130341_wp, 0.153689_wp, &
                      -0.079683_wp, -0.664796_wp, 0.150971_wp, 0.166179_wp, 0.143524_wp, &
                      -0.594149_wp, 0.74459_wp, -0.554425_wp, -0.659074_wp, -0.082647_wp, &
                      0.149176_wp, 0.126455_wp, 0.117716_wp, 0.148326_wp, 0.164576_wp, &
                      0.142436_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07658_wp, 0.072309_wp, 0.029919_wp, 0.026078_wp, 0.026329_wp, &
                      0.023964_wp, -0.078776_wp, 0.026172_wp, 0.025352_wp, 0.029547_wp, &
                      -0.208658_wp, 0.072328_wp, -0.076582_wp, -0.078772_wp, 0.023985_wp, &
                      0.026328_wp, 0.026077_wp, 0.029908_wp, 0.026173_wp, 0.025354_wp, &
                      0.029552_wp]
         case ("resp")
            charges = [ &
                      -0.570247_wp, 0.746404_wp, 0.130904_wp, 0.127277_wp, 0.149807_wp, &
                      -0.065209_wp, -0.660248_wp, 0.152441_wp, 0.162732_wp, 0.147288_wp, &
                      -0.642462_wp, 0.746873_wp, -0.554213_wp, -0.654482_wp, -0.068145_wp, &
                      0.145304_wp, 0.123423_wp, 0.125484_wp, 0.149766_wp, 0.161117_wp, &
                      0.146187_wp]
         case ("mbis")
            charges = [ &
                      -0.542768_wp, 0.31049_wp, 0.157387_wp, 0.153924_wp, 0.156702_wp, &
                      0.066153_wp, -0.555805_wp, 0.149176_wp, 0.157453_wp, 0.158178_wp, &
                      -0.421802_wp, 0.310466_wp, -0.542747_wp, -0.555796_wp, 0.066196_wp, &
                      0.15669_wp, 0.153922_wp, 0.157361_wp, 0.149178_wp, 0.15746_wp, &
                      0.158187_wp]
         case ("chelpg")
            charges = [ &
                      -0.589069_wp, 0.769505_wp, 0.134851_wp, 0.131218_wp, 0.154305_wp, &
                      -0.070085_wp, -0.683031_wp, 0.157675_wp, 0.167461_wp, 0.152386_wp, &
                      -0.650586_wp, 0.769914_wp, -0.57301_wp, -0.677129_wp, -0.07303_wp, &
                      0.149803_wp, 0.127348_wp, 0.129456_wp, 0.154982_wp, 0.16581_wp, &
                      0.151225_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.077948_wp, 0.072048_wp, 0.032534_wp, 0.025164_wp, 0.024801_wp, &
                      0.027282_wp, -0.079257_wp, 0.027678_wp, 0.024061_wp, 0.031737_wp, &
                      -0.216232_wp, 0.072067_wp, -0.077949_wp, -0.079254_wp, 0.027305_wp, &
                      0.024799_wp, 0.025164_wp, 0.032523_wp, 0.027677_wp, 0.024063_wp, &
                      0.031743_wp]
         case ("resp")
            charges = [ &
                      -0.577985_wp, 0.757807_wp, 0.137459_wp, 0.126923_wp, 0.148628_wp, &
                      -0.059026_wp, -0.666434_wp, 0.155646_wp, 0.161542_wp, 0.151903_wp, &
                      -0.673116_wp, 0.758322_wp, -0.56206_wp, -0.660592_wp, -0.061949_wp, &
                      0.144137_wp, 0.123082_wp, 0.132084_wp, 0.15295_wp, 0.159907_wp, &
                      0.15077_wp]
         case ("mbis")
            charges = [ &
                      -0.54963_wp, 0.312734_wp, 0.163413_wp, 0.153682_wp, 0.155418_wp, &
                      0.072351_wp, -0.560834_wp, 0.152454_wp, 0.155855_wp, 0.162727_wp, &
                      -0.436365_wp, 0.31271_wp, -0.549609_wp, -0.560823_wp, 0.072397_wp, &
                      0.155404_wp, 0.153682_wp, 0.163387_wp, 0.152454_wp, 0.15586_wp, &
                      0.162737_wp]
         case ("chelpg")
            charges = [ &
                      -0.596828_wp, 0.780815_wp, 0.141419_wp, 0.130872_wp, 0.153135_wp, &
                      -0.063866_wp, -0.689228_wp, 0.160882_wp, 0.166283_wp, 0.15701_wp, &
                      -0.681166_wp, 0.781274_wp, -0.580882_wp, -0.68325_wp, -0.066798_wp, &
                      0.148645_wp, 0.127017_wp, 0.13607_wp, 0.158167_wp, 0.164612_wp, &
                      0.155817_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (108)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074646_wp, -0.006355_wp, -0.039994_wp, -0.039576_wp, -0.038996_wp, &
                      -0.040732_wp, -0.006763_wp, -0.072971_wp, 0.022481_wp, 0.022449_wp, &
                      0.020262_wp, 0.017745_wp, 0.019581_wp, 0.017934_wp, 0.018608_wp, &
                      0.021265_wp, 0.020737_wp, 0.017607_wp, 0.017909_wp, 0.020404_wp, &
                      0.016969_wp, 0.02236_wp, 0.02212_wp, 0.021524_wp]
         case ("resp")
            charges = [ &
                      -0.468078_wp, 0.354108_wp, -0.22397_wp, 0.064569_wp, -0.088196_wp, &
                      -0.05496_wp, 0.305254_wp, -0.550334_wp, 0.100123_wp, 0.100482_wp, &
                      0.110062_wp, -0.072504_wp, 0.047444_wp, 0.022013_wp, -0.005418_wp, &
                      0.000613_wp, 0.014535_wp, 0.023586_wp, 0.013051_wp, 0.008964_wp, &
                      -0.059636_wp, 0.124408_wp, 0.11895_wp, 0.114933_wp]
         case ("mbis")
            charges = [ &
                      -0.441218_wp, -0.020624_wp, -0.238823_wp, -0.193603_wp, -0.217928_wp, &
                      -0.220425_wp, -0.014381_wp, -0.463466_wp, 0.129573_wp, 0.134416_wp, &
                      0.128769_wp, 0.089804_wp, 0.106402_wp, 0.101821_wp, 0.104939_wp, &
                      0.101713_wp, 0.105147_wp, 0.101815_wp, 0.105838_wp, 0.105049_wp, &
                      0.091046_wp, 0.133611_wp, 0.13287_wp, 0.137581_wp]
         case ("chelpg")
            charges = [ &
                      -0.492185_wp, 0.375857_wp, -0.245223_wp, 0.080296_wp, -0.103677_wp, &
                      -0.05076_wp, 0.311073_wp, -0.571087_wp, 0.105514_wp, 0.105675_wp, &
                      0.115966_wp, -0.076144_wp, 0.051627_wp, 0.025672_wp, -0.008431_wp, &
                      -0.001328_wp, 0.017803_wp, 0.026821_wp, 0.012007_wp, 0.00914_wp, &
                      -0.061436_wp, 0.129688_wp, 0.1237_wp, 0.119434_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.075819_wp, -0.006249_wp, -0.040569_wp, -0.040224_wp, -0.039593_wp, &
                      -0.041248_wp, -0.006598_wp, -0.074212_wp, 0.022495_wp, 0.02245_wp, &
                      0.021162_wp, 0.018171_wp, 0.019755_wp, 0.018277_wp, 0.019108_wp, &
                      0.021197_wp, 0.020709_wp, 0.01792_wp, 0.018578_wp, 0.02048_wp, &
                      0.01763_wp, 0.022395_wp, 0.022401_wp, 0.021707_wp]
         case ("resp")
            charges = [ &
                      -0.485332_wp, 0.370764_wp, -0.235119_wp, 0.069711_wp, -0.090632_wp, &
                      -0.056166_wp, 0.318202_wp, -0.571298_wp, 0.102803_wp, 0.103373_wp, &
                      0.114965_wp, -0.075914_wp, 0.049296_wp, 0.023321_wp, -0.005884_wp, &
                      -0.000997_wp, 0.013782_wp, 0.024309_wp, 0.013731_wp, 0.008283_wp, &
                      -0.061593_wp, 0.128358_wp, 0.12327_wp, 0.118767_wp]
         case ("mbis")
            charges = [ &
                      -0.448962_wp, -0.01827_wp, -0.242191_wp, -0.195592_wp, -0.219586_wp, &
                      -0.223664_wp, -0.011873_wp, -0.471373_wp, 0.130978_wp, 0.135793_wp, &
                      0.131889_wp, 0.090415_wp, 0.107422_wp, 0.103066_wp, 0.106395_wp, &
                      0.102098_wp, 0.105552_wp, 0.102863_wp, 0.107719_wp, 0.105786_wp, &
                      0.09207_wp, 0.13516_wp, 0.134839_wp, 0.139391_wp]
         case ("chelpg")
            charges = [ &
                      -0.509657_wp, 0.392778_wp, -0.256743_wp, 0.085891_wp, -0.106404_wp, &
                      -0.051833_wp, 0.323944_wp, -0.59215_wp, 0.108238_wp, 0.108607_wp, &
                      0.120933_wp, -0.079583_wp, 0.053545_wp, 0.027043_wp, -0.009_wp, &
                      -0.003028_wp, 0.017098_wp, 0.027595_wp, 0.012657_wp, 0.008454_wp, &
                      -0.063384_wp, 0.133668_wp, 0.128044_wp, 0.123287_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.076998_wp, -0.006091_wp, -0.041129_wp, -0.040873_wp, -0.040176_wp, &
                      -0.041743_wp, -0.006369_wp, -0.07546_wp, 0.022402_wp, 0.022372_wp, &
                      0.022197_wp, 0.018641_wp, 0.01986_wp, 0.01868_wp, 0.019682_wp, &
                      0.020991_wp, 0.020575_wp, 0.01826_wp, 0.019344_wp, 0.020472_wp, &
                      0.018429_wp, 0.02233_wp, 0.022677_wp, 0.02185_wp]
         case ("resp")
            charges = [ &
                      -0.503053_wp, 0.388315_wp, -0.246568_wp, 0.075288_wp, -0.092892_wp, &
                      -0.057152_wp, 0.332141_wp, -0.592978_wp, 0.105304_wp, 0.106156_wp, &
                      0.120228_wp, -0.079533_wp, 0.051007_wp, 0.024703_wp, -0.006335_wp, &
                      -0.003022_wp, 0.012717_wp, 0.024978_wp, 0.014455_wp, 0.007272_wp, &
                      -0.063557_wp, 0.132227_wp, 0.127692_wp, 0.122608_wp]
         case ("mbis")
            charges = [ &
                      -0.45673_wp, -0.015816_wp, -0.245538_wp, -0.197551_wp, -0.221097_wp, &
                      -0.226905_wp, -0.009315_wp, -0.479237_wp, 0.132184_wp, 0.137023_wp, &
                      0.135282_wp, 0.091082_wp, 0.108306_wp, 0.104417_wp, 0.107995_wp, &
                      0.102218_wp, 0.105733_wp, 0.103947_wp, 0.109779_wp, 0.106348_wp, &
                      0.093332_wp, 0.136529_wp, 0.136805_wp, 0.141133_wp]
         case ("chelpg")
            charges = [ &
                      -0.52759_wp, 0.410591_wp, -0.268557_wp, 0.091926_wp, -0.108946_wp, &
                      -0.052668_wp, 0.337801_wp, -0.613917_wp, 0.11078_wp, 0.11143_wp, &
                      0.126257_wp, -0.08323_wp, 0.05532_wp, 0.028485_wp, -0.009558_wp, &
                      -0.005146_wp, 0.016076_wp, 0.02831_wp, 0.013345_wp, 0.007433_wp, &
                      -0.065341_wp, 0.137564_wp, 0.132489_wp, 0.127146_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (109)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.062899_wp, 0.036625_wp, 0.042125_wp, 0.034099_wp, -0.049942_wp, &
                      -0.049948_wp, -0.062902_wp, 0.042127_wp, 0.034093_wp, 0.036622_wp]
         case ("resp")
            charges = [ &
                      -0.319868_wp, 0.120499_wp, 0.157114_wp, 0.130519_wp, -0.088784_wp, &
                      -0.087414_wp, -0.318601_wp, 0.156887_wp, 0.130131_wp, 0.119518_wp]
         case ("mbis")
            charges = [ &
                      -0.449793_wp, 0.17791_wp, 0.180462_wp, 0.168847_wp, -0.077417_wp, &
                      -0.077415_wp, -0.449776_wp, 0.180454_wp, 0.168826_wp, 0.177901_wp]
         case ("chelpg")
            charges = [ &
                      -0.327766_wp, 0.122653_wp, 0.159277_wp, 0.132776_wp, -0.087467_wp, &
                      -0.086092_wp, -0.326423_wp, 0.159031_wp, 0.132361_wp, 0.121651_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.058753_wp, 0.040498_wp, 0.050578_wp, 0.042346_wp, -0.07466_wp, &
                      -0.074669_wp, -0.058756_wp, 0.05058_wp, 0.042341_wp, 0.040496_wp]
         case ("resp")
            charges = [ &
                      -0.282086_wp, 0.117592_wp, 0.163167_wp, 0.13877_wp, -0.138033_wp, &
                      -0.13645_wp, -0.280925_wp, 0.162916_wp, 0.138409_wp, 0.116642_wp]
         case ("mbis")
            charges = [ &
                      -0.448967_wp, 0.183901_wp, 0.194262_wp, 0.182403_wp, -0.111587_wp, &
                      -0.111591_wp, -0.448951_wp, 0.194254_wp, 0.182384_wp, 0.183891_wp]
         case ("chelpg")
            charges = [ &
                      -0.289822_wp, 0.119705_wp, 0.165287_wp, 0.140981_wp, -0.136749_wp, &
                      -0.135162_wp, -0.288584_wp, 0.165017_wp, 0.140593_wp, 0.118734_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.057993_wp, 0.041013_wp, 0.05179_wp, 0.043701_wp, -0.078501_wp, &
                      -0.07851_wp, -0.057997_wp, 0.051792_wp, 0.043696_wp, 0.04101_wp]
         case ("resp")
            charges = [ &
                      -0.274193_wp, 0.116506_wp, 0.163426_wp, 0.139848_wp, -0.146191_wp, &
                      -0.144566_wp, -0.273084_wp, 0.16318_wp, 0.139503_wp, 0.115571_wp]
         case ("mbis")
            charges = [ &
                      -0.448206_wp, 0.184568_wp, 0.196117_wp, 0.18452_wp, -0.116987_wp, &
                      -0.116991_wp, -0.448189_wp, 0.196108_wp, 0.184501_wp, 0.184558_wp]
         case ("chelpg")
            charges = [ &
                      -0.281888_wp, 0.11861_wp, 0.165534_wp, 0.14205_wp, -0.144915_wp, &
                      -0.143286_wp, -0.280703_wp, 0.165268_wp, 0.141677_wp, 0.117653_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (110)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.339776_wp, 0.183465_wp, -0.069525_wp, -0.074633_wp, -0.007571_wp, &
                      -0.010728_wp, 0.03201_wp, 0.039034_wp, 0.031834_wp, 0.041208_wp, &
                      0.029248_wp, 0.029207_wp, 0.036495_wp, 0.043346_wp, 0.036391_wp]
         case ("resp")
            charges = [ &
                      -0.567506_wp, 0.535478_wp, -0.483169_wp, 0.068547_wp, -0.403206_wp, &
                      -0.425416_wp, 0.136449_wp, 0.13753_wp, 0.135804_wp, 0.168485_wp, &
                      0.137139_wp, 0.137848_wp, 0.132296_wp, 0.15677_wp, 0.132951_wp]
         case ("mbis")
            charges = [ &
                      -0.622355_wp, 0.716119_wp, -0.634151_wp, -0.28154_wp, -0.224022_wp, &
                      -0.224197_wp, 0.1105_wp, 0.139667_wp, 0.110059_wp, 0.160582_wp, &
                      0.104796_wp, 0.104777_wp, 0.170404_wp, 0.198986_wp, 0.170381_wp]
         case ("chelpg")
            charges = [ &
                      -0.570621_wp, 0.547963_wp, -0.507037_wp, 0.075023_wp, -0.422526_wp, &
                      -0.445449_wp, 0.141589_wp, 0.142414_wp, 0.140931_wp, 0.173298_wp, &
                      0.142409_wp, 0.143121_wp, 0.137773_wp, 0.162669_wp, 0.138445_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.404036_wp, 0.182179_wp, -0.06929_wp, -0.066556_wp, -0.001955_wp, &
                      -0.011428_wp, 0.040766_wp, 0.047823_wp, 0.04061_wp, 0.030068_wp, &
                      0.038422_wp, 0.038436_wp, 0.048833_wp, 0.037354_wp, 0.048779_wp]
         case ("resp")
            charges = [ &
                      -0.674457_wp, 0.564671_wp, -0.50059_wp, 0.077592_wp, -0.449608_wp, &
                      -0.42299_wp, 0.161228_wp, 0.162665_wp, 0.160673_wp, 0.151364_wp, &
                      0.152707_wp, 0.153487_wp, 0.156848_wp, 0.14878_wp, 0.15763_wp]
         case ("mbis")
            charges = [ &
                      -0.717613_wp, 0.75076_wp, -0.664685_wp, -0.27342_wp, -0.238059_wp, &
                      -0.239892_wp, 0.125659_wp, 0.155116_wp, 0.12525_wp, 0.14494_wp, &
                      0.123581_wp, 0.12366_wp, 0.195847_wp, 0.192948_wp, 0.195911_wp]
         case ("chelpg")
            charges = [ &
                      -0.67758_wp, 0.577178_wp, -0.524644_wp, 0.084327_wp, -0.469469_wp, &
                      -0.443276_wp, 0.166499_wp, 0.167686_wp, 0.16593_wp, 0.156224_wp, &
                      0.158038_wp, 0.158821_wp, 0.162367_wp, 0.154731_wp, 0.163168_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.4063_wp, 0.182072_wp, -0.069274_wp, -0.066238_wp, -0.001721_wp, &
                      -0.011442_wp, 0.041095_wp, 0.048148_wp, 0.04094_wp, 0.029585_wp, &
                      0.038769_wp, 0.038785_wp, 0.049288_wp, 0.037061_wp, 0.049237_wp]
         case ("resp")
            charges = [ &
                      -0.678194_wp, 0.565451_wp, -0.500968_wp, 0.077905_wp, -0.451156_wp, &
                      -0.422615_wp, 0.162121_wp, 0.163567_wp, 0.161568_wp, 0.150515_wp, &
                      0.153241_wp, 0.154025_wp, 0.157721_wp, 0.148309_wp, 0.15851_wp]
         case ("mbis")
            charges = [ &
                      -0.720952_wp, 0.751828_wp, -0.665639_wp, -0.273074_wp, -0.23848_wp, &
                      -0.240349_wp, 0.126205_wp, 0.155666_wp, 0.125797_wp, 0.144204_wp, &
                      0.124271_wp, 0.124355_wp, 0.196764_wp, 0.192575_wp, 0.196832_wp]
         case ("chelpg")
            charges = [ &
                      -0.68132_wp, 0.577955_wp, -0.525028_wp, 0.08465_wp, -0.471036_wp, &
                      -0.442908_wp, 0.167396_wp, 0.168593_wp, 0.16683_wp, 0.155377_wp, &
                      0.158575_wp, 0.159361_wp, 0.163243_wp, 0.154262_wp, 0.164049_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (111)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.333518_wp, 0.153825_wp, 0.035648_wp, -0.066225_wp, -0.003852_wp, &
                      -0.006046_wp, 0.034129_wp, 0.041816_wp, 0.034156_wp, 0.044758_wp, &
                      0.032649_wp, 0.032659_wp]
         case ("resp")
            charges = [ &
                      -0.495508_wp, 0.219916_wp, 0.074796_wp, 0.178818_wp, -0.505968_wp, &
                      -0.318282_wp, 0.161563_wp, 0.155244_wp, 0.161229_wp, 0.148492_wp, &
                      0.110028_wp, 0.109671_wp]
         case ("mbis")
            charges = [ &
                      -0.575581_wp, 0.498669_wp, 0.031113_wp, -0.237275_wp, -0.233423_wp, &
                      -0.218911_wp, 0.116327_wp, 0.126789_wp, 0.116356_wp, 0.15902_wp, &
                      0.108426_wp, 0.10849_wp]
         case ("chelpg")
            charges = [ &
                      -0.496825_wp, 0.222396_wp, 0.073638_wp, 0.188029_wp, -0.52475_wp, &
                      -0.336285_wp, 0.166352_wp, 0.15974_wp, 0.166005_wp, 0.152825_wp, &
                      0.114621_wp, 0.114256_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.401578_wp, 0.153091_wp, 0.051052_wp, -0.055982_wp, 0.001816_wp, &
                      -0.005576_wp, 0.042922_wp, 0.051159_wp, 0.042957_wp, 0.035044_wp, &
                      0.042498_wp, 0.042596_wp]
         case ("resp")
            charges = [ &
                      -0.60982_wp, 0.256273_wp, 0.097896_wp, 0.18635_wp, -0.549085_wp, &
                      -0.323543_wp, 0.185656_wp, 0.178394_wp, 0.185317_wp, 0.135212_wp, &
                      0.128769_wp, 0.128581_wp]
         case ("mbis")
            charges = [ &
                      -0.671635_wp, 0.51083_wp, 0.059706_wp, -0.221412_wp, -0.251387_wp, &
                      -0.234465_wp, 0.132043_wp, 0.143128_wp, 0.132083_wp, 0.145124_wp, &
                      0.127877_wp, 0.128107_wp]
         case ("chelpg")
            charges = [ &
                      -0.611137_wp, 0.258746_wp, 0.09674_wp, 0.195686_wp, -0.568139_wp, &
                      -0.341745_wp, 0.190515_wp, 0.182952_wp, 0.190163_wp, 0.139589_wp, &
                      0.133414_wp, 0.133218_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.403974_wp, 0.153017_wp, 0.051569_wp, -0.055584_wp, 0.002046_wp, &
                      -0.005548_wp, 0.043253_wp, 0.051482_wp, 0.043288_wp, 0.034605_wp, &
                      0.042872_wp, 0.042974_wp]
         case ("resp")
            charges = [ &
                      -0.61389_wp, 0.257678_wp, 0.098639_wp, 0.186403_wp, -0.55032_wp, &
                      -0.323498_wp, 0.186493_wp, 0.17913_wp, 0.186154_wp, 0.134491_wp, &
                      0.12945_wp, 0.12927_wp]
         case ("mbis")
            charges = [ &
                      -0.675018_wp, 0.511204_wp, 0.060667_wp, -0.220809_wp, -0.251928_wp, &
                      -0.234918_wp, 0.132609_wp, 0.143672_wp, 0.132649_wp, 0.144439_wp, &
                      0.128598_wp, 0.128834_wp]
         case ("chelpg")
            charges = [ &
                      -0.61521_wp, 0.260147_wp, 0.097484_wp, 0.195742_wp, -0.569381_wp, &
                      -0.341705_wp, 0.191354_wp, 0.18369_wp, 0.191003_wp, 0.138869_wp, &
                      0.134098_wp, 0.133909_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (112)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.063906_wp, 0.045836_wp, 0.051147_wp, 0.041052_wp, 0.260696_wp, &
                      -0.40895_wp, -0.063898_wp, 0.041035_wp, 0.045935_wp, 0.051049_wp]
         case ("resp")
            charges = [ &
                      -0.423149_wp, 0.176721_wp, 0.179453_wp, 0.14491_wp, 0.32071_wp, &
                      -0.473252_wp, -0.443882_wp, 0.151562_wp, 0.181937_wp, 0.184989_wp]
         case ("mbis")
            charges = [ &
                      -0.565551_wp, 0.197864_wp, 0.207856_wp, 0.186218_wp, 0.539063_wp, &
                      -0.591846_wp, -0.565569_wp, 0.186263_wp, 0.198067_wp, 0.207632_wp]
         case ("chelpg")
            charges = [ &
                      -0.435273_wp, 0.180209_wp, 0.182582_wp, 0.148322_wp, 0.325054_wp, &
                      -0.473398_wp, -0.455998_wp, 0.154974_wp, 0.185411_wp, 0.188118_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.054833_wp, 0.052805_wp, 0.055101_wp, 0.059976_wp, 0.260932_wp, &
                      -0.486971_wp, -0.05484_wp, 0.059975_wp, 0.052765_wp, 0.055086_wp]
         case ("resp")
            charges = [ &
                      -0.415272_wp, 0.18583_wp, 0.180361_wp, 0.176829_wp, 0.34839_wp, &
                      -0.600716_wp, -0.435068_wp, 0.183334_wp, 0.190631_wp, 0.185682_wp]
         case ("mbis")
            charges = [ &
                      -0.57136_wp, 0.208735_wp, 0.213054_wp, 0.217022_wp, 0.562902_wp, &
                      -0.697774_wp, -0.571338_wp, 0.217093_wp, 0.208701_wp, 0.212961_wp]
         case ("chelpg")
            charges = [ &
                      -0.427507_wp, 0.189347_wp, 0.183517_wp, 0.180273_wp, 0.352775_wp, &
                      -0.60086_wp, -0.447285_wp, 0.186775_wp, 0.194131_wp, 0.188835_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.054555_wp, 0.052982_wp, 0.055159_wp, 0.060572_wp, 0.260893_wp, &
                      -0.489149_wp, -0.054562_wp, 0.060572_wp, 0.052937_wp, 0.055147_wp]
         case ("resp")
            charges = [ &
                      -0.41476_wp, 0.185982_wp, 0.18022_wp, 0.177801_wp, 0.349045_wp, &
                      -0.604366_wp, -0.434518_wp, 0.184298_wp, 0.190767_wp, 0.185532_wp]
         case ("mbis")
            charges = [ &
                      -0.57139_wp, 0.208985_wp, 0.213087_wp, 0.21797_wp, 0.563475_wp, &
                      -0.700747_wp, -0.571366_wp, 0.218043_wp, 0.208942_wp, 0.212999_wp]
         case ("chelpg")
            charges = [ &
                      -0.427_wp, 0.189502_wp, 0.183379_wp, 0.181246_wp, 0.353428_wp, &
                      -0.604513_wp, -0.446739_wp, 0.18774_wp, 0.194269_wp, 0.188689_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (113)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.048402_wp, 0.084971_wp, -0.046731_wp, -0.03082_wp, -0.045585_wp, &
                      -0.031237_wp, -0.145628_wp, 0.085052_wp, -0.050323_wp, -0.031365_wp, &
                      -0.046567_wp, -0.030882_wp, -0.047666_wp, 0.042112_wp, 0.042419_wp, &
                      0.038621_wp, 0.038779_wp, 0.038448_wp, 0.037468_wp, 0.03537_wp, &
                      0.038265_wp, 0.038628_wp, 0.03507_wp]
         case ("resp")
            charges = [ &
                      -0.284764_wp, 0.457311_wp, -0.303046_wp, -0.078019_wp, -0.217674_wp, &
                      -0.085755_wp, -0.428591_wp, 0.450427_wp, -0.280666_wp, -0.093703_wp, &
                      -0.212164_wp, -0.079443_wp, -0.303294_wp, 0.170921_wp, 0.169871_wp, &
                      0.150883_wp, 0.131951_wp, 0.132423_wp, 0.148382_wp, 0.145064_wp, &
                      0.135111_wp, 0.131535_wp, 0.14324_wp]
         case ("mbis")
            charges = [ &
                      -0.302254_wp, 0.35843_wp, -0.286022_wp, -0.083339_wp, -0.188437_wp, &
                      -0.082503_wp, -0.340472_wp, 0.360656_wp, -0.306482_wp, -0.080449_wp, &
                      -0.191167_wp, -0.081924_wp, -0.287841_wp, 0.169495_wp, 0.169627_wp, &
                      0.176139_wp, 0.133972_wp, 0.133738_wp, 0.17663_wp, 0.142644_wp, &
                      0.133271_wp, 0.133632_wp, 0.142655_wp]
         case ("chelpg")
            charges = [ &
                      -0.308076_wp, 0.483616_wp, -0.325436_wp, -0.060539_wp, -0.235533_wp, &
                      -0.067149_wp, -0.439669_wp, 0.475325_wp, -0.302958_wp, -0.076179_wp, &
                      -0.228665_wp, -0.063246_wp, -0.3243_wp, 0.175487_wp, 0.174648_wp, &
                      0.155224_wp, 0.129286_wp, 0.129537_wp, 0.152582_wp, 0.148344_wp, &
                      0.132394_wp, 0.129071_wp, 0.146237_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.050771_wp, 0.079624_wp, -0.051483_wp, -0.034227_wp, -0.049193_wp, &
                      -0.034071_wp, -0.153478_wp, 0.079796_wp, -0.052678_wp, -0.034333_wp, &
                      -0.050344_wp, -0.034355_wp, -0.052574_wp, 0.045822_wp, 0.046188_wp, &
                      0.043109_wp, 0.044982_wp, 0.044823_wp, 0.041919_wp, 0.04106_wp, &
                      0.044616_wp, 0.044848_wp, 0.040717_wp]
         case ("resp")
            charges = [ &
                      -0.309566_wp, 0.472681_wp, -0.32469_wp, -0.085223_wp, -0.235767_wp, &
                      -0.090559_wp, -0.462533_wp, 0.469635_wp, -0.30783_wp, -0.09786_wp, &
                      -0.231368_wp, -0.085761_wp, -0.326706_wp, 0.185072_wp, 0.183953_wp, &
                      0.170591_wp, 0.146748_wp, 0.146421_wp, 0.169131_wp, 0.160043_wp, &
                      0.14898_wp, 0.146253_wp, 0.158354_wp]
         case ("mbis")
            charges = [ &
                      -0.309814_wp, 0.355883_wp, -0.300447_wp, -0.092314_wp, -0.201435_wp, &
                      -0.092484_wp, -0.351184_wp, 0.358561_wp, -0.313878_wp, -0.090755_wp, &
                      -0.204505_wp, -0.090808_wp, -0.302798_wp, 0.180124_wp, 0.18024_wp, &
                      0.18636_wp, 0.14774_wp, 0.14771_wp, 0.18679_wp, 0.156112_wp, &
                      0.147296_wp, 0.147454_wp, 0.156149_wp]
         case ("chelpg")
            charges = [ &
                      -0.332776_wp, 0.498914_wp, -0.346915_wp, -0.067978_wp, -0.253414_wp, &
                      -0.072126_wp, -0.473637_wp, 0.494551_wp, -0.330126_wp, -0.080395_wp, &
                      -0.247771_wp, -0.069692_wp, -0.347645_wp, 0.189627_wp, 0.188702_wp, &
                      0.174921_wp, 0.144137_wp, 0.143571_wp, 0.173339_wp, 0.163282_wp, &
                      0.146277_wp, 0.143822_wp, 0.161332_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.051829_wp, 0.076895_wp, -0.053846_wp, -0.035468_wp, -0.050297_wp, &
                      -0.034868_wp, -0.157272_wp, 0.07711_wp, -0.053714_wp, -0.035191_wp, &
                      -0.051528_wp, -0.035631_wp, -0.055004_wp, 0.047038_wp, 0.047429_wp, &
                      0.044921_wp, 0.047699_wp, 0.04773_wp, 0.043736_wp, 0.043684_wp, &
                      0.047512_wp, 0.047571_wp, 0.043318_wp]
         case ("resp")
            charges = [ &
                      -0.32187_wp, 0.478766_wp, -0.334813_wp, -0.087701_wp, -0.242553_wp, &
                      -0.091264_wp, -0.478308_wp, 0.477923_wp, -0.321546_wp, -0.098102_wp, &
                      -0.238818_wp, -0.087667_wp, -0.337836_wp, 0.190906_wp, 0.189735_wp, &
                      0.180066_wp, 0.153045_wp, 0.152326_wp, 0.17926_wp, 0.166398_wp, &
                      0.154791_wp, 0.152471_wp, 0.164789_wp]
         case ("mbis")
            charges = [ &
                      -0.313033_wp, 0.353988_wp, -0.307113_wp, -0.095574_wp, -0.206331_wp, &
                      -0.09593_wp, -0.356082_wp, 0.356869_wp, -0.317_wp, -0.094356_wp, &
                      -0.209557_wp, -0.094043_wp, -0.309694_wp, 0.184178_wp, 0.184287_wp, &
                      0.190536_wp, 0.153574_wp, 0.153733_wp, 0.19097_wp, 0.161939_wp, &
                      0.153343_wp, 0.15331_wp, 0.161983_wp]
         case ("chelpg")
            charges = [ &
                      -0.345107_wp, 0.50502_wp, -0.357024_wp, -0.070501_wp, -0.260168_wp, &
                      -0.072836_wp, -0.48944_wp, 0.502904_wp, -0.343917_wp, -0.08059_wp, &
                      -0.255241_wp, -0.071593_wp, -0.358807_wp, 0.195466_wp, 0.194481_wp, &
                      0.184407_wp, 0.150448_wp, 0.149479_wp, 0.183488_wp, 0.169632_wp, &
                      0.152081_wp, 0.150045_wp, 0.167772_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (114)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.182401_wp, 0.007927_wp, 0.007926_wp, 0.007405_wp, 0.023192_wp, &
                      0.088645_wp, 0.023194_wp, 0.007405_wp, -0.035192_wp, 0.023322_wp, &
                      0.019598_wp, -0.069037_wp, -0.035186_wp, 0.019606_wp, 0.023329_wp, &
                      -0.069044_wp, 0.025089_wp, 0.023086_wp, 0.021488_wp, 0.025088_wp, &
                      0.021486_wp, 0.023082_wp]
         case ("resp")
            charges = [ &
                      -0.891668_wp, 0.358504_wp, 0.357153_wp, -0.062719_wp, -0.005446_wp, &
                      0.3159_wp, -0.005404_wp, -0.062032_wp, 0.147606_wp, -0.005024_wp, &
                      -0.009368_wp, -0.419241_wp, 0.144258_wp, -0.008983_wp, -0.004849_wp, &
                      -0.410964_wp, 0.110335_wp, 0.08665_wp, 0.086952_wp, 0.10837_wp, &
                      0.085053_wp, 0.084918_wp]
         case ("mbis")
            charges = [ &
                      -0.60507_wp, 0.032481_wp, 0.03245_wp, 0.04558_wp, 0.086508_wp, &
                      0.316541_wp, 0.086497_wp, 0.045603_wp, -0.217725_wp, 0.120853_wp, &
                      0.102649_wp, -0.412473_wp, -0.217677_wp, 0.102667_wp, 0.120824_wp, &
                      -0.412507_wp, 0.132154_wp, 0.128357_wp, 0.125888_wp, 0.132159_wp, &
                      0.125887_wp, 0.128363_wp]
         case ("chelpg")
            charges = [ &
                      -0.902263_wp, 0.370447_wp, 0.369577_wp, -0.066017_wp, -0.008481_wp, &
                      0.316304_wp, -0.008564_wp, -0.065441_wp, 0.155507_wp, -0.006306_wp, &
                      -0.010976_wp, -0.438018_wp, 0.151838_wp, -0.010551_wp, -0.006089_wp, &
                      -0.42945_wp, 0.115119_wp, 0.090764_wp, 0.091254_wp, 0.113101_wp, &
                      0.08928_wp, 0.088967_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.193633_wp, 0.006805_wp, 0.006805_wp, 0.010326_wp, 0.023882_wp, &
                      0.094853_wp, 0.023885_wp, 0.010326_wp, -0.0364_wp, 0.020356_wp, &
                      0.022337_wp, -0.07015_wp, -0.0364_wp, 0.022341_wp, 0.020368_wp, &
                      -0.070155_wp, 0.025789_wp, 0.023801_wp, 0.022643_wp, 0.025788_wp, &
                      0.022642_wp, 0.023798_wp]
         case ("resp")
            charges = [ &
                      -0.987522_wp, 0.396144_wp, 0.395187_wp, -0.062183_wp, -0.009793_wp, &
                      0.350396_wp, -0.009837_wp, -0.061515_wp, 0.152556_wp, -0.012347_wp, &
                      -0.005944_wp, -0.432284_wp, 0.148986_wp, -0.005571_wp, -0.012134_wp, &
                      -0.423803_wp, 0.114664_wp, 0.088839_wp, 0.089207_wp, 0.112676_wp, &
                      0.087244_wp, 0.087036_wp]
         case ("mbis")
            charges = [ &
                      -0.652324_wp, 0.039075_wp, 0.039054_wp, 0.050737_wp, 0.087977_wp, &
                      0.338665_wp, 0.087966_wp, 0.050759_wp, -0.221439_wp, 0.116792_wp, &
                      0.108934_wp, -0.420049_wp, -0.221421_wp, 0.108949_wp, 0.116777_wp, &
                      -0.420068_wp, 0.134741_wp, 0.130931_wp, 0.129137_wp, 0.134744_wp, &
                      0.129136_wp, 0.130936_wp]
         case ("chelpg")
            charges = [ &
                      -0.998194_wp, 0.408208_wp, 0.407741_wp, -0.06552_wp, -0.012865_wp, &
                      0.350786_wp, -0.013036_wp, -0.064966_wp, 0.160563_wp, -0.013659_wp, &
                      -0.007576_wp, -0.451215_wp, 0.156669_wp, -0.007164_wp, -0.013405_wp, &
                      -0.442441_wp, 0.119486_wp, 0.092987_wp, 0.093538_wp, 0.117445_wp, &
                      0.0915_wp, 0.091118_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.200944_wp, 0.006178_wp, 0.006178_wp, 0.012919_wp, 0.023884_wp, &
                      0.098714_wp, 0.023888_wp, 0.012919_wp, -0.037183_wp, 0.017737_wp, &
                      0.024295_wp, -0.070653_wp, -0.037185_wp, 0.024293_wp, 0.017753_wp, &
                      -0.070656_wp, 0.026149_wp, 0.024076_wp, 0.023713_wp, 0.026147_wp, &
                      0.023713_wp, 0.024074_wp]
         case ("resp")
            charges = [ &
                      -1.051043_wp, 0.421356_wp, 0.420611_wp, -0.060377_wp, -0.013578_wp, &
                      0.372857_wp, -0.013668_wp, -0.059709_wp, 0.154353_wp, -0.018165_wp, &
                      -0.002899_wp, -0.438786_wp, 0.150665_wp, -0.002543_wp, -0.017935_wp, &
                      -0.430157_wp, 0.11697_wp, 0.089526_wp, 0.090948_wp, 0.114961_wp, &
                      0.088939_wp, 0.087674_wp]
         case ("mbis")
            charges = [ &
                      -0.683381_wp, 0.043649_wp, 0.043635_wp, 0.055352_wp, 0.088024_wp, &
                      0.352615_wp, 0.088014_wp, 0.055372_wp, -0.223527_wp, 0.112835_wp, &
                      0.11333_wp, -0.424008_wp, -0.223528_wp, 0.113342_wp, 0.11283_wp, &
                      -0.424016_wp, 0.136047_wp, 0.132062_wp, 0.131624_wp, 0.136047_wp, &
                      0.131623_wp, 0.132067_wp]
         case ("chelpg")
            charges = [ &
                      -1.061745_wp, 0.43349_wp, 0.433239_wp, -0.063735_wp, -0.016673_wp, &
                      0.373232_wp, -0.016893_wp, -0.063183_wp, 0.162419_wp, -0.019497_wp, &
                      -0.004544_wp, -0.457794_wp, 0.158406_wp, -0.004148_wp, -0.019225_wp, &
                      -0.448871_wp, 0.12181_wp, 0.093691_wp, 0.095293_wp, 0.119748_wp, &
                      0.093208_wp, 0.091773_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (115)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070365_wp, -0.034214_wp, -0.036429_wp, -0.036262_wp, -0.036209_wp, &
                      -0.036241_wp, -0.036241_wp, -0.036207_wp, -0.036261_wp, -0.036431_wp, &
                      -0.034205_wp, -0.070371_wp, 0.021442_wp, 0.021419_wp, 0.023249_wp, &
                      0.019065_wp, 0.019081_wp, 0.018181_wp, 0.018211_wp, 0.018204_wp, &
                      0.018209_wp, 0.01817_wp, 0.018188_wp, 0.018171_wp, 0.018162_wp, &
                      0.018162_wp, 0.01817_wp, 0.018189_wp, 0.018169_wp, 0.018204_wp, &
                      0.018207_wp, 0.018211_wp, 0.018175_wp, 0.019078_wp, 0.019067_wp, &
                      0.021426_wp, 0.021432_wp, 0.023258_wp]
         case ("resp")
            charges = [ &
                      -0.345895_wp, 0.210916_wp, 0.000326_wp, -0.053506_wp, 0.096339_wp, &
                      0.004168_wp, -0.005815_wp, 0.100405_wp, -0.05693_wp, 0.005147_wp, &
                      0.228772_wp, -0.355344_wp, 0.073944_wp, 0.074027_wp, 0.080228_wp, &
                      -0.033955_wp, -0.034014_wp, -0.00142_wp, -0.001459_wp, -0.000625_wp, &
                      -0.000933_wp, -0.025394_wp, -0.024724_wp, -0.007934_wp, -0.007883_wp, &
                      -0.005705_wp, -0.006158_wp, -0.025185_wp, -0.02668_wp, -0.000167_wp, &
                      -4.3e-05_wp, -0.005374_wp, -0.002844_wp, -0.040541_wp, -0.038966_wp, &
                      0.075029_wp, 0.076028_wp, 0.082168_wp]
         case ("mbis")
            charges = [ &
                      -0.421884_wp, -0.160064_wp, -0.211283_wp, -0.213974_wp, -0.208426_wp, &
                      -0.209916_wp, -0.209888_wp, -0.20844_wp, -0.21393_wp, -0.211293_wp, &
                      -0.160051_wp, -0.421824_wp, 0.127566_wp, 0.127544_wp, 0.129692_wp, &
                      0.099735_wp, 0.09986_wp, 0.105526_wp, 0.105612_wp, 0.105239_wp, &
                      0.105315_wp, 0.104766_wp, 0.104793_wp, 0.104987_wp, 0.104935_wp, &
                      0.10496_wp, 0.104956_wp, 0.104815_wp, 0.104741_wp, 0.1053_wp, &
                      0.105235_wp, 0.105611_wp, 0.105511_wp, 0.099884_wp, 0.099698_wp, &
                      0.127591_wp, 0.127473_wp, 0.129688_wp]
         case ("chelpg")
            charges = [ &
                      -0.365599_wp, 0.226128_wp, -0.000573_wp, -0.067877_wp, 0.116563_wp, &
                      -0.00338_wp, -0.014927_wp, 0.121384_wp, -0.071848_wp, 0.00472_wp, &
                      0.244446_wp, -0.375385_wp, 0.078334_wp, 0.078421_wp, 0.084658_wp, &
                      -0.036592_wp, -0.036652_wp, -0.001497_wp, -0.001526_wp, 0.001809_wp, &
                      0.001523_wp, -0.029044_wp, -0.028368_wp, -0.006941_wp, -0.006904_wp, &
                      -0.004408_wp, -0.00484_wp, -0.028926_wp, -0.030451_wp, 0.002398_wp, &
                      0.002495_wp, -0.005629_wp, -0.003036_wp, -0.043332_wp, -0.041759_wp, &
                      0.079484_wp, 0.080446_wp, 0.086684_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071728_wp, -0.034697_wp, -0.036906_wp, -0.036775_wp, -0.03675_wp, &
                      -0.036788_wp, -0.036787_wp, -0.036748_wp, -0.036772_wp, -0.036908_wp, &
                      -0.034688_wp, -0.071734_wp, 0.021714_wp, 0.021697_wp, 0.023491_wp, &
                      0.019426_wp, 0.019437_wp, 0.01853_wp, 0.018555_wp, 0.018514_wp, &
                      0.018515_wp, 0.018447_wp, 0.018466_wp, 0.018443_wp, 0.01844_wp, &
                      0.018436_wp, 0.018448_wp, 0.018462_wp, 0.01845_wp, 0.018508_wp, &
                      0.018519_wp, 0.018548_wp, 0.018529_wp, 0.01943_wp, 0.019432_wp, &
                      0.021693_wp, 0.021715_wp, 0.023499_wp]
         case ("resp")
            charges = [ &
                      -0.357018_wp, 0.221076_wp, 0.00234_wp, -0.052327_wp, 0.102606_wp, &
                      0.007311_wp, -0.002424_wp, 0.106585_wp, -0.055836_wp, 0.007544_wp, &
                      0.23916_wp, -0.366455_wp, 0.075922_wp, 0.076022_wp, 0.082313_wp, &
                      -0.036129_wp, -0.036219_wp, -0.002485_wp, -0.002573_wp, -0.001987_wp, &
                      -0.002332_wp, -0.027782_wp, -0.027092_wp, -0.009798_wp, -0.009696_wp, &
                      -0.007574_wp, -0.008028_wp, -0.02755_wp, -0.029056_wp, -0.001573_wp, &
                      -0.001427_wp, -0.006639_wp, -0.004007_wp, -0.042873_wp, -0.041221_wp, &
                      0.076958_wp, 0.078013_wp, 0.084251_wp]
         case ("mbis")
            charges = [ &
                      -0.429227_wp, -0.16037_wp, -0.213303_wp, -0.216014_wp, -0.210376_wp, &
                      -0.211876_wp, -0.211847_wp, -0.210388_wp, -0.215968_wp, -0.21331_wp, &
                      -0.160356_wp, -0.42917_wp, 0.129446_wp, 0.129434_wp, 0.131528_wp, &
                      0.10066_wp, 0.100773_wp, 0.106626_wp, 0.106702_wp, 0.106273_wp, &
                      0.106342_wp, 0.105745_wp, 0.105774_wp, 0.105963_wp, 0.105922_wp, &
                      0.105939_wp, 0.105941_wp, 0.105786_wp, 0.105727_wp, 0.106323_wp, &
                      0.106274_wp, 0.106686_wp, 0.106623_wp, 0.10079_wp, 0.100631_wp, &
                      0.129463_wp, 0.129375_wp, 0.131522_wp]
         case ("chelpg")
            charges = [ &
                      -0.376855_wp, 0.236363_wp, 0.001567_wp, -0.066789_wp, 0.123036_wp, &
                      -0.000189_wp, -0.011387_wp, 0.127739_wp, -0.070884_wp, 0.007328_wp, &
                      0.25487_wp, -0.38662_wp, 0.080341_wp, 0.080444_wp, 0.086778_wp, &
                      -0.038786_wp, -0.038877_wp, -0.0026_wp, -0.002679_wp, 0.000444_wp, &
                      0.000122_wp, -0.031498_wp, -0.0308_wp, -0.008851_wp, -0.008762_wp, &
                      -0.00634_wp, -0.006773_wp, -0.031349_wp, -0.032887_wp, 0.000993_wp, &
                      0.001114_wp, -0.00695_wp, -0.004255_wp, -0.045678_wp, -0.044028_wp, &
                      0.081439_wp, 0.082457_wp, 0.088803_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073156_wp, -0.035166_wp, -0.037386_wp, -0.037296_wp, -0.037301_wp, &
                      -0.037343_wp, -0.037342_wp, -0.037299_wp, -0.037293_wp, -0.037387_wp, &
                      -0.035157_wp, -0.073163_wp, 0.02197_wp, 0.02196_wp, 0.023688_wp, &
                      0.019818_wp, 0.019822_wp, 0.018901_wp, 0.018919_wp, 0.018846_wp, &
                      0.018842_wp, 0.018727_wp, 0.018747_wp, 0.018717_wp, 0.018722_wp, &
                      0.018714_wp, 0.018729_wp, 0.018736_wp, 0.018735_wp, 0.018834_wp, &
                      0.018854_wp, 0.018902_wp, 0.018907_wp, 0.019811_wp, 0.019829_wp, &
                      0.021942_wp, 0.021985_wp, 0.023695_wp]
         case ("resp")
            charges = [ &
                      -0.368702_wp, 0.232552_wp, 0.004919_wp, -0.050532_wp, 0.10991_wp, &
                      0.011261_wp, 0.001801_wp, 0.113791_wp, -0.054128_wp, 0.010555_wp, &
                      0.250888_wp, -0.378124_wp, 0.077898_wp, 0.078018_wp, 0.084368_wp, &
                      -0.038647_wp, -0.038774_wp, -0.003839_wp, -0.003986_wp, -0.003689_wp, &
                      -0.004076_wp, -0.030634_wp, -0.029922_wp, -0.012097_wp, -0.011936_wp, &
                      -0.009876_wp, -0.010332_wp, -0.030378_wp, -0.031892_wp, -0.003328_wp, &
                      -0.003154_wp, -0.008225_wp, -0.005469_wp, -0.04557_wp, -0.043831_wp, &
                      0.078877_wp, 0.08_wp, 0.086303_wp]
         case ("mbis")
            charges = [ &
                      -0.436852_wp, -0.160596_wp, -0.215357_wp, -0.218102_wp, -0.212357_wp, &
                      -0.213864_wp, -0.213834_wp, -0.212368_wp, -0.218054_wp, -0.21536_wp, &
                      -0.160581_wp, -0.436797_wp, 0.131361_wp, 0.131363_wp, 0.133346_wp, &
                      0.101629_wp, 0.101728_wp, 0.107767_wp, 0.107828_wp, 0.107349_wp, &
                      0.107408_wp, 0.106735_wp, 0.106765_wp, 0.106948_wp, 0.106921_wp, &
                      0.10693_wp, 0.106938_wp, 0.106764_wp, 0.106727_wp, 0.107386_wp, &
                      0.107355_wp, 0.107793_wp, 0.107777_wp, 0.101737_wp, 0.101609_wp, &
                      0.131366_wp, 0.131316_wp, 0.133339_wp]
         case ("chelpg")
            charges = [ &
                      -0.38866_wp, 0.24789_wp, 0.004313_wp, -0.06502_wp, 0.1305_wp, &
                      0.00386_wp, -0.006944_wp, 0.135072_wp, -0.06925_wp, 0.010604_wp, &
                      0.266607_wp, -0.398403_wp, 0.082341_wp, 0.082465_wp, 0.088867_wp, &
                      -0.041324_wp, -0.041451_wp, -0.004006_wp, -0.004144_wp, -0.001279_wp, &
                      -0.001642_wp, -0.034413_wp, -0.033694_wp, -0.011214_wp, -0.011066_wp, &
                      -0.008726_wp, -0.009162_wp, -0.034236_wp, -0.035782_wp, -0.000777_wp, &
                      -0.000629_wp, -0.008609_wp, -0.00579_wp, -0.048388_wp, -0.04665_wp, &
                      0.083381_wp, 0.084466_wp, 0.090891_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (116)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069786_wp, -0.020187_wp, -0.106279_wp, 0.028033_wp, 0.026649_wp, &
                      0.028337_wp, 0.037725_wp, 0.031315_wp, 0.044191_wp]
         case ("resp")
            charges = [ &
                      -0.251916_wp, 0.100482_wp, -0.346721_wp, 0.065464_wp, 0.079667_wp, &
                      0.070036_wp, 0.099255_wp, 0.008944_wp, 0.174789_wp]
         case ("mbis")
            charges = [ &
                      -0.405934_wp, -0.196399_wp, -0.27841_wp, 0.132081_wp, 0.139096_wp, &
                      0.149081_wp, 0.156779_wp, 0.141664_wp, 0.16204_wp]
         case ("chelpg")
            charges = [ &
                      -0.263998_wp, 0.107198_wp, -0.347813_wp, 0.068413_wp, 0.082513_wp, &
                      0.072976_wp, 0.098112_wp, 0.007605_wp, 0.174995_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06986_wp, -0.01661_wp, -0.139505_wp, 0.03281_wp, 0.02947_wp, &
                      0.026928_wp, 0.042977_wp, 0.037128_wp, 0.056662_wp]
         case ("resp")
            charges = [ &
                      -0.27017_wp, 0.153558_wp, -0.428785_wp, 0.073495_wp, 0.086591_wp, &
                      0.069072_wp, 0.098459_wp, 0.00812_wp, 0.209659_wp]
         case ("mbis")
            charges = [ &
                      -0.417823_wp, -0.180905_wp, -0.344167_wp, 0.141548_wp, 0.145473_wp, &
                      0.1481_wp, 0.164293_wp, 0.149673_wp, 0.193807_wp]
         case ("chelpg")
            charges = [ &
                      -0.282847_wp, 0.160954_wp, -0.429974_wp, 0.076577_wp, 0.089552_wp, &
                      0.072126_wp, 0.097145_wp, 0.006605_wp, 0.209861_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069748_wp, -0.0157_wp, -0.14709_wp, 0.033982_wp, 0.030161_wp, &
                      0.026439_wp, 0.043965_wp, 0.038561_wp, 0.059429_wp]
         case ("resp")
            charges = [ &
                      -0.27438_wp, 0.166364_wp, -0.447807_wp, 0.075504_wp, 0.088321_wp, &
                      0.068563_wp, 0.097655_wp, 0.008198_wp, 0.217581_wp]
         case ("mbis")
            charges = [ &
                      -0.419938_wp, -0.177142_wp, -0.359053_wp, 0.14374_wp, 0.146899_wp, &
                      0.147448_wp, 0.165585_wp, 0.151636_wp, 0.200825_wp]
         case ("chelpg")
            charges = [ &
                      -0.287156_wp, 0.173866_wp, -0.449008_wp, 0.078608_wp, 0.091301_wp, &
                      0.071636_wp, 0.096315_wp, 0.006656_wp, 0.217781_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (117)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074704_wp, 0.023072_wp, 0.044527_wp, 0.028964_wp, 0.029941_wp, &
                      0.033314_wp, 0.019806_wp, -0.25205_wp, 0.147131_wp]
         case ("resp")
            charges = [ &
                      -0.378852_wp, 0.093202_wp, 0.355537_wp, 0.109091_wp, 0.081714_wp, &
                      0.036056_wp, -0.036222_wp, -0.643946_wp, 0.383419_wp]
         case ("mbis")
            charges = [ &
                      -0.526687_wp, 0.137761_wp, 0.142189_wp, 0.157968_wp, 0.152461_wp, &
                      0.089317_wp, 0.052265_wp, -0.620022_wp, 0.414748_wp]
         case ("chelpg")
            charges = [ &
                      -0.392571_wp, 0.096267_wp, 0.365628_wp, 0.112322_wp, 0.084943_wp, &
                      0.034141_wp, -0.038413_wp, -0.646164_wp, 0.383847_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074822_wp, 0.029647_wp, 0.043573_wp, 0.026255_wp, 0.034553_wp, &
                      0.032575_wp, 0.026267_wp, -0.281078_wp, 0.163031_wp]
         case ("resp")
            charges = [ &
                      -0.384774_wp, 0.104_wp, 0.403277_wp, 0.102469_wp, 0.086981_wp, &
                      0.02492_wp, -0.031878_wp, -0.743091_wp, 0.438096_wp]
         case ("mbis")
            charges = [ &
                      -0.538301_wp, 0.150744_wp, 0.154275_wp, 0.154234_wp, 0.16152_wp, &
                      0.087794_wp, 0.06289_wp, -0.694719_wp, 0.461562_wp]
         case ("chelpg")
            charges = [ &
                      -0.398623_wp, 0.10709_wp, 0.41351_wp, 0.105727_wp, 0.09024_wp, &
                      0.022967_wp, -0.034103_wp, -0.745326_wp, 0.438518_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.074805_wp, 0.030028_wp, 0.043531_wp, 0.026056_wp, 0.034805_wp, &
                      0.032459_wp, 0.026637_wp, -0.282502_wp, 0.16379_wp]
         case ("resp")
            charges = [ &
                      -0.384879_wp, 0.104618_wp, 0.405522_wp, 0.101966_wp, 0.087263_wp, &
                      0.024231_wp, -0.031506_wp, -0.747978_wp, 0.440762_wp]
         case ("mbis")
            charges = [ &
                      -0.538762_wp, 0.151466_wp, 0.154893_wp, 0.153905_wp, 0.161986_wp, &
                      0.087572_wp, 0.063507_wp, -0.698356_wp, 0.463789_wp]
         case ("chelpg")
            charges = [ &
                      -0.398732_wp, 0.10771_wp, 0.415757_wp, 0.105224_wp, 0.090523_wp, &
                      0.022276_wp, -0.03373_wp, -0.750214_wp, 0.441186_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (118)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.229937_wp, -0.303462_wp, -0.149469_wp, -0.06602_wp, 0.04576_wp, &
                      0.047881_wp, 0.040525_wp, 0.055434_wp, -0.067632_wp, 0.034487_wp, &
                      0.03482_wp, 0.031427_wp, 0.034745_wp, 0.031564_wp]
         case ("resp")
            charges = [ &
                      0.85275_wp, -0.593894_wp, -0.45679_wp, -0.657681_wp, 0.187601_wp, &
                      0.184723_wp, 0.17834_wp, 0.341998_wp, -0.430975_wp, 0.017202_wp, &
                      0.017964_wp, 0.127784_wp, 0.103191_wp, 0.127787_wp]
         case ("mbis")
            charges = [ &
                      0.840494_wp, -0.602161_wp, -0.414763_wp, -0.665133_wp, 0.197991_wp, &
                      0.192701_wp, 0.189735_wp, 0.111202_wp, -0.492594_wp, 0.092041_wp, &
                      0.092313_wp, 0.152984_wp, 0.151932_wp, 0.153257_wp]
         case ("chelpg")
            charges = [ &
                      0.860659_wp, -0.595814_wp, -0.460986_wp, -0.674336_wp, 0.191622_wp, &
                      0.188645_wp, 0.182802_wp, 0.355401_wp, -0.446348_wp, 0.014163_wp, &
                      0.014944_wp, 0.131426_wp, 0.106467_wp, 0.131356_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.234801_wp, -0.338092_wp, -0.149884_wp, -0.065359_wp, 0.04856_wp, &
                      0.055816_wp, 0.046731_wp, 0.055292_wp, -0.067662_wp, 0.036697_wp, &
                      0.036882_wp, 0.032784_wp, 0.040576_wp, 0.032856_wp]
         case ("resp")
            charges = [ &
                      0.904039_wp, -0.66377_wp, -0.487624_wp, -0.681321_wp, 0.19583_wp, &
                      0.200711_wp, 0.193312_wp, 0.372932_wp, -0.442752_wp, 0.016733_wp, &
                      0.016958_wp, 0.130525_wp, 0.114124_wp, 0.130302_wp]
         case ("mbis")
            charges = [ &
                      0.878198_wp, -0.658151_wp, -0.423458_wp, -0.691842_wp, 0.205783_wp, &
                      0.209213_wp, 0.203348_wp, 0.110625_wp, -0.502256_wp, 0.096405_wp, &
                      0.096379_wp, 0.156223_wp, 0.163147_wp, 0.156385_wp]
         case ("chelpg")
            charges = [ &
                      0.911967_wp, -0.665692_wp, -0.491835_wp, -0.69805_wp, 0.19987_wp, &
                      0.204651_wp, 0.197793_wp, 0.386431_wp, -0.458242_wp, 0.013673_wp, &
                      0.013916_wp, 0.134194_wp, 0.117427_wp, 0.133898_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.235792_wp, -0.346499_wp, -0.149662_wp, -0.065101_wp, 0.049055_wp, &
                      0.057752_wp, 0.04838_wp, 0.055255_wp, -0.067536_wp, 0.037051_wp, &
                      0.037182_wp, 0.033125_wp, 0.042036_wp, 0.033168_wp]
         case ("resp")
            charges = [ &
                      0.914904_wp, -0.680502_wp, -0.49381_wp, -0.685767_wp, 0.197249_wp, &
                      0.204463_wp, 0.196988_wp, 0.379695_wp, -0.444624_wp, 0.016473_wp, &
                      0.016495_wp, 0.131026_wp, 0.11671_wp, 0.130701_wp]
         case ("mbis")
            charges = [ &
                      0.886547_wp, -0.671619_wp, -0.424867_wp, -0.697624_wp, 0.207239_wp, &
                      0.213117_wp, 0.206771_wp, 0.110603_wp, -0.504057_wp, 0.097125_wp, &
                      0.096992_wp, 0.156912_wp, 0.165839_wp, 0.15702_wp]
         case ("chelpg")
            charges = [ &
                      0.922841_wp, -0.682437_wp, -0.498024_wp, -0.702517_wp, 0.201295_wp, &
                      0.20841_wp, 0.201474_wp, 0.39322_wp, -0.460141_wp, 0.013408_wp, &
                      0.013449_wp, 0.1347_wp, 0.12002_wp, 0.134302_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (119)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.20764_wp, -0.288903_wp, -0.14419_wp, 0.055378_wp, 0.058651_wp, &
                      -0.066161_wp, 0.037358_wp, 0.037359_wp, 0.036266_wp, 0.033301_wp, &
                      0.0333_wp]
         case ("resp")
            charges = [ &
                      0.587724_wp, -0.528918_wp, -0.42965_wp, 0.028906_wp, 0.374992_wp, &
                      -0.496037_wp, 0.026741_wp, 0.027053_wp, 0.115077_wp, 0.147085_wp, &
                      0.147029_wp]
         case ("mbis")
            charges = [ &
                      0.604444_wp, -0.543004_wp, -0.388604_wp, 0.049723_wp, 0.117568_wp, &
                      -0.496201_wp, 0.094928_wp, 0.09493_wp, 0.154281_wp, 0.155971_wp, &
                      0.155962_wp]
         case ("chelpg")
            charges = [ &
                      0.589785_wp, -0.530085_wp, -0.4326_wp, 0.028607_wp, 0.387442_wp, &
                      -0.512857_wp, 0.024249_wp, 0.024584_wp, 0.118778_wp, 0.151093_wp, &
                      0.151005_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.215464_wp, -0.325911_wp, -0.144282_wp, 0.069736_wp, 0.059095_wp, &
                      -0.066168_wp, 0.040451_wp, 0.04045_wp, 0.042457_wp, 0.034353_wp, &
                      0.034354_wp]
         case ("resp")
            charges = [ &
                      0.640769_wp, -0.602719_wp, -0.457819_wp, 0.044158_wp, 0.400312_wp, &
                      -0.504638_wp, 0.028349_wp, 0.028672_wp, 0.126672_wp, 0.148145_wp, &
                      0.148097_wp]
         case ("mbis")
            charges = [ &
                      0.630819_wp, -0.600848_wp, -0.396905_wp, 0.071388_wp, 0.116676_wp, &
                      -0.506102_wp, 0.100672_wp, 0.100668_wp, 0.166202_wp, 0.158717_wp, &
                      0.158712_wp]
         case ("chelpg")
            charges = [ &
                      0.642767_wp, -0.603867_wp, -0.460748_wp, 0.043882_wp, 0.412782_wp, &
                      -0.521548_wp, 0.025859_wp, 0.026206_wp, 0.130394_wp, 0.152176_wp, &
                      0.152097_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.216697_wp, -0.331983_wp, -0.144036_wp, 0.072048_wp, 0.059165_wp, &
                      -0.066084_wp, 0.040831_wp, 0.040829_wp, 0.043506_wp, 0.034513_wp, &
                      0.034514_wp]
         case ("resp")
            charges = [ &
                      0.649235_wp, -0.614841_wp, -0.461704_wp, 0.046622_wp, 0.404032_wp, &
                      -0.50512_wp, 0.028419_wp, 0.028743_wp, 0.12848_wp, 0.14809_wp, &
                      0.148044_wp]
         case ("mbis")
            charges = [ &
                      0.634995_wp, -0.610337_wp, -0.397779_wp, 0.074862_wp, 0.116572_wp, &
                      -0.507332_wp, 0.10138_wp, 0.101375_wp, 0.168142_wp, 0.159063_wp, &
                      0.159059_wp]
         case ("chelpg")
            charges = [ &
                      0.65123_wp, -0.615993_wp, -0.464634_wp, 0.046348_wp, 0.41651_wp, &
                      -0.522043_wp, 0.025929_wp, 0.026277_wp, 0.132205_wp, 0.152124_wp, &
                      0.152046_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (120)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.082653_wp, -0.157564_wp, 0.053543_wp, -0.064718_wp, -0.053447_wp, &
                      -0.066965_wp, -0.0324_wp, -0.034546_wp, -0.054133_wp, 0.028061_wp, &
                      0.02804_wp, 0.030956_wp, 0.038441_wp, 0.032623_wp, 0.032209_wp, &
                      0.032626_wp, 0.036854_wp, 0.035444_wp, 0.032323_wp]
         case ("resp")
            charges = [ &
                      0.413244_wp, -0.415007_wp, 0.394668_wp, -0.256862_wp, -0.32391_wp, &
                      -0.471675_wp, -0.051845_wp, -0.109944_wp, -0.228172_wp, -0.006287_wp, &
                      -0.00719_wp, 0.110723_wp, 0.168328_wp, 0.137241_wp, 0.106129_wp, &
                      0.136653_wp, 0.123927_wp, 0.138137_wp, 0.141841_wp]
         case ("mbis")
            charges = [ &
                      0.384977_wp, -0.3451_wp, 0.120682_wp, -0.355969_wp, -0.298503_wp, &
                      -0.497459_wp, -0.07314_wp, -0.058795_wp, -0.211817_wp, 0.073979_wp, &
                      0.074112_wp, 0.160334_wp, 0.16594_wp, 0.156343_wp, 0.148592_wp, &
                      0.156379_wp, 0.130225_wp, 0.127288_wp, 0.141929_wp]
         case ("chelpg")
            charges = [ &
                      0.438423_wp, -0.424018_wp, 0.408867_wp, -0.278936_wp, -0.345481_wp, &
                      -0.486023_wp, -0.035425_wp, -0.093437_wp, -0.243768_wp, -0.009548_wp, &
                      -0.010446_wp, 0.115231_wp, 0.173176_wp, 0.140751_wp, 0.108991_wp, &
                      0.140171_wp, 0.121274_wp, 0.135582_wp, 0.144618_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.076011_wp, -0.167993_wp, 0.054414_wp, -0.066488_wp, -0.061654_wp, &
                      -0.068568_wp, -0.038086_wp, -0.039175_wp, -0.060425_wp, 0.033995_wp, &
                      0.03403_wp, 0.039645_wp, 0.040936_wp, 0.031897_wp, 0.03776_wp, &
                      0.031901_wp, 0.042687_wp, 0.041795_wp, 0.037315_wp]
         case ("resp")
            charges = [ &
                      0.427325_wp, -0.455503_wp, 0.414243_wp, -0.276585_wp, -0.348369_wp, &
                      -0.491697_wp, -0.061678_wp, -0.117803_wp, -0.254113_wp, 0.002247_wp, &
                      0.001527_wp, 0.131395_wp, 0.180532_wp, 0.138828_wp, 0.120201_wp, &
                      0.138297_wp, 0.139718_wp, 0.153185_wp, 0.158251_wp]
         case ("mbis")
            charges = [ &
                      0.383571_wp, -0.362055_wp, 0.117499_wp, -0.365752_wp, -0.320839_wp, &
                      -0.510068_wp, -0.083138_wp, -0.072074_wp, -0.229603_wp, 0.08498_wp, &
                      0.085222_wp, 0.17848_wp, 0.176704_wp, 0.157087_wp, 0.160455_wp, &
                      0.157124_wp, 0.144358_wp, 0.142334_wp, 0.155712_wp]
         case ("chelpg")
            charges = [ &
                      0.452259_wp, -0.464418_wp, 0.428237_wp, -0.298439_wp, -0.369577_wp, &
                      -0.506112_wp, -0.045707_wp, -0.101553_wp, -0.26935_wp, -0.000943_wp, &
                      -0.001657_wp, 0.135861_wp, 0.185298_wp, 0.142365_wp, 0.12309_wp, &
                      0.141843_wp, 0.137171_wp, 0.150677_wp, 0.160955_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.073288_wp, -0.172186_wp, 0.054806_wp, -0.066709_wp, -0.065192_wp, &
                      -0.069115_wp, -0.040313_wp, -0.040438_wp, -0.062592_wp, 0.036399_wp, &
                      0.036459_wp, 0.04325_wp, 0.041348_wp, 0.031289_wp, 0.040189_wp, &
                      0.031292_wp, 0.044702_wp, 0.044348_wp, 0.039173_wp]
         case ("resp")
            charges = [ &
                      0.432253_wp, -0.471741_wp, 0.421851_wp, -0.283549_wp, -0.358119_wp, &
                      -0.498477_wp, -0.065422_wp, -0.1195_wp, -0.263638_wp, 0.005752_wp, &
                      0.005108_wp, 0.139714_wp, 0.184467_wp, 0.138625_wp, 0.126134_wp, &
                      0.138116_wp, 0.145425_wp, 0.158734_wp, 0.164268_wp]
         case ("mbis")
            charges = [ &
                      0.382516_wp, -0.368836_wp, 0.116184_wp, -0.368501_wp, -0.329894_wp, &
                      -0.514113_wp, -0.086729_wp, -0.076244_wp, -0.235982_wp, 0.089399_wp, &
                      0.089688_wp, 0.185604_wp, 0.180055_wp, 0.156664_wp, 0.165435_wp, &
                      0.1567_wp, 0.149352_wp, 0.147977_wp, 0.160723_wp]
         case ("chelpg")
            charges = [ &
                      0.457104_wp, -0.480622_wp, 0.435735_wp, -0.305348_wp, -0.379191_wp, &
                      -0.512903_wp, -0.049603_wp, -0.103346_wp, -0.278751_wp, 0.002597_wp, &
                      0.001959_wp, 0.144176_wp, 0.189209_wp, 0.14217_wp, 0.129035_wp, &
                      0.14167_wp, 0.142915_wp, 0.156243_wp, 0.166953_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (121)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.041703_wp, -0.037017_wp, -0.041342_wp, 0.008179_wp, -0.041348_wp, &
                      -0.037018_wp, 0.034298_wp, 0.034841_wp, 0.031996_wp, -0.03323_wp, &
                      0.03199_wp, 0.034841_wp, 0.026063_wp, -0.068351_wp, 0.026063_wp, &
                      0.024494_wp, 0.023618_wp, 0.023621_wp]
         case ("resp")
            charges = [ &
                      -0.09496_wp, -0.168372_wp, -0.158997_wp, 0.043255_wp, -0.161665_wp, &
                      -0.167067_wp, 0.117578_wp, 0.133988_wp, 0.129255_wp, 0.084946_wp, &
                      0.130417_wp, 0.133783_wp, 0.001847_wp, -0.13634_wp, 0.002148_wp, &
                      0.029357_wp, 0.040414_wp, 0.040412_wp]
         case ("mbis")
            charges = [ &
                      -0.154941_wp, -0.121205_wp, -0.198083_wp, 0.112738_wp, -0.198103_wp, &
                      -0.121194_wp, 0.13596_wp, 0.132567_wp, 0.136254_wp, -0.214145_wp, &
                      0.13624_wp, 0.132565_wp, 0.116986_wp, -0.410834_wp, 0.116991_wp, &
                      0.128767_wp, 0.13471_wp, 0.134721_wp]
         case ("chelpg")
            charges = [ &
                      -0.100726_wp, -0.161669_wp, -0.167872_wp, 0.051058_wp, -0.170553_wp, &
                      -0.160369_wp, 0.118465_wp, 0.132986_wp, 0.131032_wp, 0.09092_wp, &
                      0.132199_wp, 0.132784_wp, 0.000213_wp, -0.145625_wp, 0.000515_wp, &
                      0.031615_wp, 0.042515_wp, 0.042513_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.046918_wp, -0.041891_wp, -0.0452_wp, 0.005371_wp, -0.045208_wp, &
                      -0.041893_wp, 0.037904_wp, 0.038667_wp, 0.036379_wp, -0.033442_wp, &
                      0.036373_wp, 0.038666_wp, 0.028409_wp, -0.069279_wp, 0.028406_wp, &
                      0.026345_wp, 0.023655_wp, 0.023652_wp]
         case ("resp")
            charges = [ &
                      -0.103003_wp, -0.184403_wp, -0.168204_wp, 0.032097_wp, -0.170876_wp, &
                      -0.183089_wp, 0.12744_wp, 0.145486_wp, 0.141647_wp, 0.09301_wp, &
                      0.142789_wp, 0.145284_wp, 0.003161_wp, -0.13168_wp, 0.003458_wp, &
                      0.030505_wp, 0.038189_wp, 0.03819_wp]
         case ("mbis")
            charges = [ &
                      -0.167255_wp, -0.132504_wp, -0.21036_wp, 0.112784_wp, -0.210385_wp, &
                      -0.132493_wp, 0.146103_wp, 0.142865_wp, 0.147393_wp, -0.217972_wp, &
                      0.147378_wp, 0.142863_wp, 0.121825_wp, -0.416916_wp, 0.121824_wp, &
                      0.133124_wp, 0.135861_wp, 0.135861_wp]
         case ("chelpg")
            charges = [ &
                      -0.108401_wp, -0.178127_wp, -0.176473_wp, 0.039092_wp, -0.179157_wp, &
                      -0.176817_wp, 0.128272_wp, 0.14455_wp, 0.143303_wp, 0.09951_wp, &
                      0.144451_wp, 0.144352_wp, 0.001432_wp, -0.141123_wp, 0.00173_wp, &
                      0.032766_wp, 0.040318_wp, 0.040319_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.051299_wp, -0.045901_wp, -0.048266_wp, 0.00307_wp, -0.048276_wp, &
                      -0.045903_wp, 0.040656_wp, 0.041666_wp, 0.040005_wp, -0.033508_wp, &
                      0.039998_wp, 0.041665_wp, 0.03044_wp, -0.069859_wp, 0.030435_wp, &
                      0.028147_wp, 0.023467_wp, 0.023458_wp]
         case ("resp")
            charges = [ &
                      -0.10982_wp, -0.197373_wp, -0.175171_wp, 0.021592_wp, -0.177836_wp, &
                      -0.196056_wp, 0.135237_wp, 0.154638_wp, 0.151766_wp, 0.100418_wp, &
                      0.152889_wp, 0.154439_wp, 0.004343_wp, -0.127206_wp, 0.004637_wp, &
                      0.031993_wp, 0.035755_wp, 0.035753_wp]
         case ("mbis")
            charges = [ &
                      -0.177448_wp, -0.141622_wp, -0.220044_wp, 0.112393_wp, -0.220074_wp, &
                      -0.14161_wp, 0.154075_wp, 0.151018_wp, 0.156467_wp, -0.220978_wp, &
                      0.156451_wp, 0.151015_wp, 0.125912_wp, -0.420933_wp, 0.125908_wp, &
                      0.137034_wp, 0.136221_wp, 0.136211_wp]
         case ("chelpg")
            charges = [ &
                      -0.11481_wp, -0.191561_wp, -0.18279_wp, 0.027737_wp, -0.185466_wp, &
                      -0.190248_wp, 0.136008_wp, 0.153776_wp, 0.153296_wp, 0.107418_wp, &
                      0.154426_wp, 0.15358_wp, 0.00253_wp, -0.136793_wp, 0.002825_wp, &
                      0.034256_wp, 0.03791_wp, 0.037908_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (122)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.026537_wp, 0.041193_wp, -0.046888_wp, 0.047529_wp, 0.099166_wp, &
                      -0.125369_wp, -0.046891_wp, 0.047527_wp, -0.026539_wp, 0.041194_wp, &
                      -0.041972_wp, 0.037575_wp]
         case ("resp")
            charges = [ &
                      -0.070565_wp, 0.132368_wp, -0.29044_wp, 0.175485_wp, 0.409674_wp, &
                      -0.240031_wp, -0.288165_wp, 0.174666_wp, -0.073341_wp, 0.133573_wp, &
                      -0.206778_wp, 0.143554_wp]
         case ("mbis")
            charges = [ &
                      -0.072637_wp, 0.135587_wp, -0.30171_wp, 0.177776_wp, 0.383893_wp, &
                      -0.218668_wp, -0.301724_wp, 0.17778_wp, -0.072667_wp, 0.135594_wp, &
                      -0.188445_wp, 0.145207_wp]
         case ("chelpg")
            charges = [ &
                      -0.052753_wp, 0.129484_wp, -0.311585_wp, 0.179704_wp, 0.431785_wp, &
                      -0.243845_wp, -0.309233_wp, 0.178866_wp, -0.055632_wp, 0.130725_wp, &
                      -0.224206_wp, 0.146691_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.029143_wp, 0.048881_wp, -0.050601_wp, 0.053753_wp, 0.090847_wp, &
                      -0.136422_wp, -0.050606_wp, 0.053751_wp, -0.029145_wp, 0.048883_wp, &
                      -0.0447_wp, 0.04449_wp]
         case ("resp")
            charges = [ &
                      -0.076658_wp, 0.149045_wp, -0.306671_wp, 0.191674_wp, 0.408449_wp, &
                      -0.257068_wp, -0.304509_wp, 0.190891_wp, -0.079328_wp, 0.150213_wp, &
                      -0.227325_wp, 0.161287_wp]
         case ("mbis")
            charges = [ &
                      -0.0821_wp, 0.151767_wp, -0.314165_wp, 0.192133_wp, 0.374993_wp, &
                      -0.229843_wp, -0.314185_wp, 0.192136_wp, -0.08213_wp, 0.151775_wp, &
                      -0.200946_wp, 0.160551_wp]
         case ("chelpg")
            charges = [ &
                      -0.059143_wp, 0.146222_wp, -0.327578_wp, 0.195851_wp, 0.430381_wp, &
                      -0.26087_wp, -0.325343_wp, 0.195051_wp, -0.061908_wp, 0.147423_wp, &
                      -0.244458_wp, 0.164372_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.029611_wp, 0.050997_wp, -0.051727_wp, 0.055223_wp, 0.088313_wp, &
                      -0.139427_wp, -0.051734_wp, 0.055221_wp, -0.029613_wp, 0.050999_wp, &
                      -0.045112_wp, 0.046456_wp]
         case ("resp")
            charges = [ &
                      -0.077739_wp, 0.153451_wp, -0.311014_wp, 0.195702_wp, 0.407023_wp, &
                      -0.26145_wp, -0.308884_wp, 0.19493_wp, -0.080378_wp, 0.154607_wp, &
                      -0.232275_wp, 0.166027_wp]
         case ("mbis")
            charges = [ &
                      -0.084132_wp, 0.156042_wp, -0.317583_wp, 0.195665_wp, 0.371887_wp, &
                      -0.232739_wp, -0.317606_wp, 0.195668_wp, -0.084162_wp, 0.15605_wp, &
                      -0.20376_wp, 0.164656_wp]
         case ("chelpg")
            charges = [ &
                      -0.060259_wp, 0.150638_wp, -0.331904_wp, 0.199876_wp, 0.42895_wp, &
                      -0.265252_wp, -0.329704_wp, 0.199086_wp, -0.062991_wp, 0.151827_wp, &
                      -0.249374_wp, 0.169107_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (123)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.331068_wp, 0.163452_wp, 0.040753_wp, -0.146821_wp, 0.133982_wp, &
                      0.139702_wp]
         case ("resp")
            charges = [ &
                      -0.527937_wp, 0.545676_wp, 0.02892_wp, -0.82138_wp, 0.36686_wp, &
                      0.407862_wp]
         case ("mbis")
            charges = [ &
                      -0.582498_wp, 0.65216_wp, 0.01998_wp, -0.913666_wp, 0.403314_wp, &
                      0.420709_wp]
         case ("chelpg")
            charges = [ &
                      -0.529779_wp, 0.552356_wp, 0.027181_wp, -0.828751_wp, 0.369039_wp, &
                      0.409954_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.407545_wp, 0.163891_wp, 0.057544_wp, -0.128358_wp, 0.160686_wp, &
                      0.153781_wp]
         case ("resp")
            charges = [ &
                      -0.648965_wp, 0.571422_wp, 0.055734_wp, -0.819677_wp, 0.408784_wp, &
                      0.432702_wp]
         case ("mbis")
            charges = [ &
                      -0.689134_wp, 0.665884_wp, 0.050105_wp, -0.911849_wp, 0.443666_wp, &
                      0.441327_wp]
         case ("chelpg")
            charges = [ &
                      -0.650758_wp, 0.57795_wp, 0.054047_wp, -0.826947_wp, 0.410934_wp, &
                      0.434774_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.40844_wp, 0.163882_wp, 0.057731_wp, -0.128105_wp, 0.161001_wp, &
                      0.153931_wp]
         case ("resp")
            charges = [ &
                      -0.650397_wp, 0.571701_wp, 0.05603_wp, -0.819516_wp, 0.409248_wp, &
                      0.432934_wp]
         case ("mbis")
            charges = [ &
                      -0.690378_wp, 0.665997_wp, 0.050444_wp, -0.911683_wp, 0.444108_wp, &
                      0.441512_wp]
         case ("chelpg")
            charges = [ &
                      -0.652198_wp, 0.578229_wp, 0.054344_wp, -0.82678_wp, 0.411397_wp, &
                      0.435007_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (124)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.17794_wp, 0.216985_wp, -0.288474_wp, 0.186524_wp, 0.062904_wp]
         case ("resp")
            charges = [ &
                      -0.541372_wp, 0.517749_wp, -0.491754_wp, 0.431873_wp, 0.083505_wp]
         case ("mbis")
            charges = [ &
                      -0.617737_wp, 0.64433_wp, -0.550463_wp, 0.464906_wp, 0.058963_wp]
         case ("chelpg")
            charges = [ &
                      -0.542363_wp, 0.519839_wp, -0.492618_wp, 0.432181_wp, 0.082961_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.177535_wp, 0.226331_wp, -0.332328_wp, 0.205487_wp, 0.078044_wp]
         case ("resp")
            charges = [ &
                      -0.57771_wp, 0.575441_wp, -0.572607_wp, 0.477513_wp, 0.097362_wp]
         case ("mbis")
            charges = [ &
                      -0.645428_wp, 0.679862_wp, -0.618594_wp, 0.503559_wp, 0.080599_wp]
         case ("chelpg")
            charges = [ &
                      -0.578689_wp, 0.577512_wp, -0.573463_wp, 0.477813_wp, 0.096828_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.177449_wp, 0.226535_wp, -0.333365_wp, 0.205892_wp, 0.078387_wp]
         case ("resp")
            charges = [ &
                      -0.578368_wp, 0.576752_wp, -0.574525_wp, 0.478473_wp, 0.097667_wp]
         case ("mbis")
            charges = [ &
                      -0.64589_wp, 0.680652_wp, -0.620199_wp, 0.504349_wp, 0.081088_wp]
         case ("chelpg")
            charges = [ &
                      -0.579363_wp, 0.578843_wp, -0.575399_wp, 0.478785_wp, 0.097135_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (125)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070355_wp, -0.034184_wp, -0.036381_wp, -0.036271_wp, -0.036381_wp, &
                      -0.034185_wp, -0.070355_wp, 0.021444_wp, 0.021442_wp, 0.023273_wp, &
                      0.019102_wp, 0.019094_wp, 0.018223_wp, 0.018226_wp, 0.01826_wp, &
                      0.018256_wp, 0.018223_wp, 0.018226_wp, 0.019097_wp, 0.019099_wp, &
                      0.021443_wp, 0.021442_wp, 0.023272_wp]
         case ("resp")
            charges = [ &
                      -0.37043_wp, 0.200592_wp, 0.072883_wp, -0.15089_wp, 0.067679_wp, &
                      0.197225_wp, -0.363813_wp, 0.079402_wp, 0.079272_wp, 0.087803_wp, &
                      -0.031154_wp, -0.030432_wp, -0.01414_wp, -0.013422_wp, 0.016134_wp, &
                      0.015136_wp, -0.011815_wp, -0.012034_wp, -0.030158_wp, -0.029537_wp, &
                      0.077492_wp, 0.077944_wp, 0.086264_wp]
         case ("mbis")
            charges = [ &
                      -0.421728_wp, -0.160127_wp, -0.209901_wp, -0.218166_wp, -0.209902_wp, &
                      -0.160128_wp, -0.421729_wp, 0.127493_wp, 0.127489_wp, 0.129688_wp, &
                      0.09987_wp, 0.099837_wp, 0.105427_wp, 0.105414_wp, 0.105629_wp, &
                      0.105628_wp, 0.10541_wp, 0.105431_wp, 0.099845_wp, 0.099861_wp, &
                      0.127491_wp, 0.12749_wp, 0.129687_wp]
         case ("chelpg")
            charges = [ &
                      -0.388853_wp, 0.209463_wp, 0.083521_wp, -0.169363_wp, 0.077975_wp, &
                      0.206257_wp, -0.38244_wp, 0.083389_wp, 0.08326_wp, 0.092538_wp, &
                      -0.032772_wp, -0.032069_wp, -0.016288_wp, -0.015593_wp, 0.019174_wp, &
                      0.018191_wp, -0.013902_wp, -0.014086_wp, -0.031795_wp, -0.031163_wp, &
                      0.081541_wp, 0.081981_wp, 0.091036_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071644_wp, -0.034635_wp, -0.036821_wp, -0.036718_wp, -0.036821_wp, &
                      -0.034635_wp, -0.071645_wp, 0.021698_wp, 0.021696_wp, 0.0235_wp, &
                      0.019442_wp, 0.019438_wp, 0.018551_wp, 0.018556_wp, 0.018585_wp, &
                      0.018582_wp, 0.018553_wp, 0.018554_wp, 0.01944_wp, 0.01944_wp, &
                      0.021698_wp, 0.021696_wp, 0.0235_wp]
         case ("resp")
            charges = [ &
                      -0.38141_wp, 0.209832_wp, 0.07734_wp, -0.153117_wp, 0.071995_wp, &
                      0.206296_wp, -0.374647_wp, 0.081334_wp, 0.081196_wp, 0.089929_wp, &
                      -0.033149_wp, -0.03242_wp, -0.015627_wp, -0.01489_wp, 0.015509_wp, &
                      0.0145_wp, -0.013242_wp, -0.013447_wp, -0.032097_wp, -0.031483_wp, &
                      0.079391_wp, 0.079843_wp, 0.088363_wp]
         case ("mbis")
            charges = [ &
                      -0.428681_wp, -0.160404_wp, -0.211787_wp, -0.220163_wp, -0.211787_wp, &
                      -0.160405_wp, -0.428681_wp, 0.129269_wp, 0.129265_wp, 0.131425_wp, &
                      0.100739_wp, 0.100712_wp, 0.106454_wp, 0.106444_wp, 0.106651_wp, &
                      0.10665_wp, 0.106441_wp, 0.106458_wp, 0.100719_wp, 0.100732_wp, &
                      0.129268_wp, 0.129266_wp, 0.131424_wp]
         case ("chelpg")
            charges = [ &
                      -0.399987_wp, 0.218747_wp, 0.088286_wp, -0.171909_wp, 0.082604_wp, &
                      0.215369_wp, -0.393432_wp, 0.085351_wp, 0.085214_wp, 0.094708_wp, &
                      -0.034784_wp, -0.034074_wp, -0.017849_wp, -0.017136_wp, 0.018587_wp, &
                      0.017594_wp, -0.015405_wp, -0.015574_wp, -0.033749_wp, -0.033124_wp, &
                      0.083471_wp, 0.083911_wp, 0.09318_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073142_wp, -0.035114_wp, -0.037303_wp, -0.037208_wp, -0.037302_wp, &
                      -0.035114_wp, -0.073143_wp, 0.021958_wp, 0.021955_wp, 0.023698_wp, &
                      0.019851_wp, 0.019851_wp, 0.018932_wp, 0.018939_wp, 0.018987_wp, &
                      0.018983_wp, 0.018936_wp, 0.018935_wp, 0.019851_wp, 0.01985_wp, &
                      0.021958_wp, 0.021955_wp, 0.023697_wp]
         case ("resp")
            charges = [ &
                      -0.394146_wp, 0.22137_wp, 0.083032_wp, -0.155246_wp, 0.077515_wp, &
                      0.217625_wp, -0.387208_wp, 0.08346_wp, 0.083311_wp, 0.092245_wp, &
                      -0.035707_wp, -0.034967_wp, -0.017601_wp, -0.016841_wp, 0.014578_wp, &
                      0.013558_wp, -0.015143_wp, -0.015331_wp, -0.034581_wp, -0.033977_wp, &
                      0.081477_wp, 0.08193_wp, 0.090647_wp]
         case ("mbis")
            charges = [ &
                      -0.436656_wp, -0.160619_wp, -0.213895_wp, -0.222417_wp, -0.213895_wp, &
                      -0.16062_wp, -0.436657_wp, 0.131262_wp, 0.131257_wp, 0.133316_wp, &
                      0.101744_wp, 0.101724_wp, 0.107621_wp, 0.107613_wp, 0.107849_wp, &
                      0.107848_wp, 0.107612_wp, 0.107623_wp, 0.101729_wp, 0.101738_wp, &
                      0.131261_wp, 0.131257_wp, 0.133315_wp]
         case ("chelpg")
            charges = [ &
                      -0.412898_wp, 0.230331_wp, 0.09435_wp, -0.174373_wp, 0.088501_wp, &
                      0.22674_wp, -0.40617_wp, 0.08751_wp, 0.087363_wp, 0.097073_wp, &
                      -0.03736_wp, -0.03664_wp, -0.019916_wp, -0.019181_wp, 0.017689_wp, &
                      0.016684_wp, -0.0174_wp, -0.017552_wp, -0.036252_wp, -0.035637_wp, &
                      0.085592_wp, 0.086032_wp, 0.095515_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (126)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070579_wp, -0.034208_wp, -0.036501_wp, -0.036348_wp, -0.036297_wp, &
                      -0.036335_wp, -0.036339_wp, -0.036341_wp, -0.036339_wp, -0.036341_wp, &
                      0.023243_wp, 0.021473_wp, 0.021475_wp, 0.019094_wp, 0.019099_wp, &
                      0.018231_wp, 0.018221_wp, 0.018237_wp, 0.018252_wp, 0.018232_wp, &
                      0.018215_wp, 0.018193_wp, 0.01821_wp, 0.018211_wp, 0.018193_wp, &
                      0.018187_wp, 0.018204_wp, 0.018198_wp, 0.018191_wp, -0.036333_wp, &
                      0.018201_wp, 0.018201_wp, -0.036298_wp, 0.018199_wp, 0.018203_wp, &
                      -0.036346_wp, 0.018224_wp, 0.018221_wp, -0.036502_wp, 0.018244_wp, &
                      0.018244_wp, -0.034207_wp, 0.018224_wp, 0.018224_wp, -0.070578_wp, &
                      0.019095_wp, 0.019098_wp, 0.023243_wp, 0.021473_wp, 0.021477_wp]
         case ("resp")
            charges = [ &
                      -0.336161_wp, 0.211918_wp, 0.010389_wp, -0.059021_wp, 0.085099_wp, &
                      0.037403_wp, -0.015489_wp, 0.052244_wp, 0.038336_wp, -0.01368_wp, &
                      0.076928_wp, 0.071306_wp, 0.071224_wp, -0.035699_wp, -0.035449_wp, &
                      -0.004325_wp, -0.004422_wp, -0.000739_wp, -0.000327_wp, -0.023177_wp, &
                      -0.022756_wp, -0.014272_wp, -0.013372_wp, -0.00647_wp, -0.008089_wp, &
                      -0.016738_wp, -0.01778_wp, -0.014004_wp, -0.013594_wp, 0.026246_wp, &
                      -0.006462_wp, -0.004851_wp, 0.094784_wp, -0.011423_wp, -0.010967_wp, &
                      -0.067373_wp, -0.025165_wp, -0.024648_wp, 0.01457_wp, 0.001268_wp, &
                      0.001087_wp, 0.213236_wp, -0.005044_wp, -0.004959_wp, -0.343702_wp, &
                      -0.036015_wp, -0.034987_wp, 0.078661_wp, 0.073097_wp, 0.073364_wp]
         case ("mbis")
            charges = [ &
                      -0.421504_wp, -0.160763_wp, -0.211792_wp, -0.214455_wp, -0.209044_wp, &
                      -0.210402_wp, -0.210458_wp, -0.210343_wp, -0.21036_wp, -0.210461_wp, &
                      0.12959_wp, 0.127551_wp, 0.127548_wp, 0.099981_wp, 0.099991_wp, &
                      0.105791_wp, 0.105789_wp, 0.105524_wp, 0.105554_wp, 0.105087_wp, &
                      0.105058_wp, 0.105211_wp, 0.10525_wp, 0.105241_wp, 0.10521_wp, &
                      0.105199_wp, 0.105211_wp, 0.105222_wp, 0.105194_wp, -0.210405_wp, &
                      0.105219_wp, 0.105232_wp, -0.209039_wp, 0.105229_wp, 0.105234_wp, &
                      -0.214456_wp, 0.105066_wp, 0.105075_wp, -0.211789_wp, 0.105533_wp, &
                      0.105545_wp, -0.160765_wp, 0.105781_wp, 0.105792_wp, -0.421496_wp, &
                      0.099988_wp, 0.099984_wp, 0.129587_wp, 0.127544_wp, 0.127555_wp]
         case ("chelpg")
            charges = [ &
                      -0.355407_wp, 0.224911_wp, 0.013075_wp, -0.073209_wp, 0.098379_wp, &
                      0.039757_wp, -0.025996_wp, 0.059984_wp, 0.044048_wp, -0.023056_wp, &
                      0.081481_wp, 0.075536_wp, 0.075486_wp, -0.037995_wp, -0.037716_wp, &
                      -0.004994_wp, -0.005099_wp, 0.001574_wp, 0.001949_wp, -0.025668_wp, &
                      -0.025266_wp, -0.014953_wp, -0.014033_wp, -0.005043_wp, -0.006644_wp, &
                      -0.018245_wp, -0.019288_wp, -0.0151_wp, -0.014695_wp, 0.027599_wp, &
                      -0.005171_wp, -0.003579_wp, 0.108549_wp, -0.011843_wp, -0.011392_wp, &
                      -0.082185_wp, -0.027732_wp, -0.027239_wp, 0.017769_wp, 0.003715_wp, &
                      0.003521_wp, 0.226199_wp, -0.005832_wp, -0.005731_wp, -0.363019_wp, &
                      -0.038327_wp, -0.037283_wp, 0.083245_wp, 0.077345_wp, 0.077618_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071954_wp, -0.034689_wp, -0.036984_wp, -0.036869_wp, -0.036841_wp, &
                      -0.036885_wp, -0.036895_wp, -0.036899_wp, -0.036897_wp, -0.036897_wp, &
                      0.023486_wp, 0.021747_wp, 0.021751_wp, 0.019458_wp, 0.019463_wp, &
                      0.01858_wp, 0.018573_wp, 0.018553_wp, 0.018564_wp, 0.01851_wp, &
                      0.018498_wp, 0.018474_wp, 0.018487_wp, 0.018484_wp, 0.01847_wp, &
                      0.018463_wp, 0.018478_wp, 0.018472_wp, 0.018468_wp, -0.036883_wp, &
                      0.018478_wp, 0.018474_wp, -0.036842_wp, 0.018478_wp, 0.018482_wp, &
                      -0.036868_wp, 0.018504_wp, 0.018502_wp, -0.036985_wp, 0.018558_wp, &
                      0.018558_wp, -0.034687_wp, 0.018573_wp, 0.018576_wp, -0.071953_wp, &
                      0.019457_wp, 0.019463_wp, 0.023485_wp, 0.021747_wp, 0.021753_wp]
         case ("resp")
            charges = [ &
                      -0.347354_wp, 0.221945_wp, 0.012764_wp, -0.058472_wp, 0.091366_wp, &
                      0.041718_wp, -0.013057_wp, 0.057228_wp, 0.042761_wp, -0.010734_wp, &
                      0.079056_wp, 0.073298_wp, 0.07322_wp, -0.037841_wp, -0.037584_wp, &
                      -0.005447_wp, -0.005555_wp, -0.00201_wp, -0.001604_wp, -0.025556_wp, &
                      -0.025146_wp, -0.016301_wp, -0.015394_wp, -0.008223_wp, -0.009857_wp, &
                      -0.018903_wp, -0.019959_wp, -0.016081_wp, -0.015671_wp, 0.030226_wp, &
                      -0.00832_wp, -0.0067_wp, 0.101053_wp, -0.013392_wp, -0.012923_wp, &
                      -0.066613_wp, -0.027548_wp, -0.027051_wp, 0.017334_wp, -0.000111_wp, &
                      -0.000281_wp, 0.223523_wp, -0.006336_wp, -0.00624_wp, -0.355174_wp, &
                      -0.038268_wp, -0.037205_wp, 0.080863_wp, 0.075139_wp, 0.075418_wp]
         case ("mbis")
            charges = [ &
                      -0.428919_wp, -0.161062_wp, -0.21383_wp, -0.216524_wp, -0.211007_wp, &
                      -0.212373_wp, -0.212442_wp, -0.212326_wp, -0.212344_wp, -0.212444_wp, &
                      0.131443_wp, 0.12945_wp, 0.129451_wp, 0.10091_wp, 0.100919_wp, &
                      0.106893_wp, 0.106897_wp, 0.106572_wp, 0.106597_wp, 0.106071_wp, &
                      0.10605_wp, 0.106207_wp, 0.106237_wp, 0.106225_wp, 0.106202_wp, &
                      0.106188_wp, 0.106198_wp, 0.106207_wp, 0.106185_wp, -0.212377_wp, &
                      0.106209_wp, 0.106217_wp, -0.211001_wp, 0.106219_wp, 0.106226_wp, &
                      -0.216526_wp, 0.106054_wp, 0.106063_wp, -0.213827_wp, 0.106578_wp, &
                      0.106591_wp, -0.161065_wp, 0.106884_wp, 0.1069_wp, -0.428911_wp, &
                      0.100915_wp, 0.100915_wp, 0.13144_wp, 0.129444_wp, 0.129458_wp]
         case ("chelpg")
            charges = [ &
                      -0.366769_wp, 0.234971_wp, 0.015739_wp, -0.072916_wp, 0.10489_wp, &
                      0.044326_wp, -0.02373_wp, 0.065223_wp, 0.048691_wp, -0.020235_wp, &
                      0.083657_wp, 0.077563_wp, 0.077517_wp, -0.040148_wp, -0.039863_wp, &
                      -0.006186_wp, -0.006303_wp, 0.000325_wp, 0.000693_wp, -0.028119_wp, &
                      -0.027728_wp, -0.017061_wp, -0.016134_wp, -0.006804_wp, -0.008421_wp, &
                      -0.020489_wp, -0.021546_wp, -0.017251_wp, -0.016847_wp, 0.031853_wp, &
                      -0.007048_wp, -0.005447_wp, 0.114999_wp, -0.013894_wp, -0.01343_wp, &
                      -0.081652_wp, -0.030173_wp, -0.029701_wp, 0.020838_wp, 0.002354_wp, &
                      0.00217_wp, 0.23651_wp, -0.007199_wp, -0.007087_wp, -0.374647_wp, &
                      -0.040594_wp, -0.039514_wp, 0.085494_wp, 0.079419_wp, 0.079704_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073349_wp, -0.035141_wp, -0.037441_wp, -0.037365_wp, -0.037372_wp, &
                      -0.037423_wp, -0.03744_wp, -0.037446_wp, -0.037443_wp, -0.037442_wp, &
                      0.023677_wp, 0.021995_wp, 0.022003_wp, 0.019839_wp, 0.019843_wp, &
                      0.018933_wp, 0.018931_wp, 0.018877_wp, 0.018884_wp, 0.018779_wp, &
                      0.018774_wp, 0.01875_wp, 0.018757_wp, 0.018745_wp, 0.018738_wp, &
                      0.01873_wp, 0.018742_wp, 0.018735_wp, 0.018736_wp, -0.037421_wp, &
                      0.018744_wp, 0.018737_wp, -0.037373_wp, 0.018751_wp, 0.018755_wp, &
                      -0.037363_wp, 0.018775_wp, 0.018774_wp, -0.037442_wp, 0.01888_wp, &
                      0.01888_wp, -0.03514_wp, 0.018928_wp, 0.018933_wp, -0.073347_wp, &
                      0.019837_wp, 0.019845_wp, 0.023677_wp, 0.021996_wp, 0.022005_wp]
         case ("resp")
            charges = [ &
                      -0.358721_wp, 0.232894_wp, 0.015671_wp, -0.057364_wp, 0.098441_wp, &
                      0.04679_wp, -0.009981_wp, 0.063042_wp, 0.04796_wp, -0.007105_wp, &
                      0.081088_wp, 0.075221_wp, 0.075149_wp, -0.040243_wp, -0.03998_wp, &
                      -0.006831_wp, -0.006949_wp, -0.003569_wp, -0.003173_wp, -0.028313_wp, &
                      -0.02791_wp, -0.018699_wp, -0.01779_wp, -0.010341_wp, -0.011988_wp, &
                      -0.021461_wp, -0.022534_wp, -0.018543_wp, -0.018131_wp, 0.034936_wp, &
                      -0.010548_wp, -0.008919_wp, 0.108123_wp, -0.015725_wp, -0.015244_wp, &
                      -0.065281_wp, -0.030306_wp, -0.029828_wp, 0.020659_wp, -0.001785_wp, &
                      -0.001945_wp, 0.234748_wp, -0.007901_wp, -0.007792_wp, -0.366845_wp, &
                      -0.04079_wp, -0.039688_wp, 0.082975_wp, 0.077118_wp, 0.077409_wp]
         case ("mbis")
            charges = [ &
                      -0.436367_wp, -0.161276_wp, -0.215814_wp, -0.21854_wp, -0.212931_wp, &
                      -0.214303_wp, -0.214384_wp, -0.214266_wp, -0.214284_wp, -0.214385_wp, &
                      0.133219_wp, 0.131321_wp, 0.131328_wp, 0.101852_wp, 0.101859_wp, &
                      0.107987_wp, 0.108_wp, 0.107619_wp, 0.107635_wp, 0.107027_wp, &
                      0.10702_wp, 0.107183_wp, 0.1072_wp, 0.107178_wp, 0.107169_wp, &
                      0.107153_wp, 0.107156_wp, 0.107165_wp, 0.107151_wp, -0.214307_wp, &
                      0.107173_wp, 0.107173_wp, -0.212925_wp, 0.10719_wp, 0.107197_wp, &
                      -0.218543_wp, 0.107015_wp, 0.107026_wp, -0.215811_wp, 0.107621_wp, &
                      0.107634_wp, -0.16128_wp, 0.107979_wp, 0.108_wp, -0.436359_wp, &
                      0.101853_wp, 0.101858_wp, 0.133217_wp, 0.131316_wp, 0.131334_wp]
         case ("chelpg")
            charges = [ &
                      -0.378287_wp, 0.245935_wp, 0.018959_wp, -0.072035_wp, 0.112184_wp, &
                      0.049665_wp, -0.020774_wp, 0.071277_wp, 0.054111_wp, -0.016693_wp, &
                      0.085735_wp, 0.079516_wp, 0.079478_wp, -0.04256_wp, -0.042268_wp, &
                      -0.007647_wp, -0.007774_wp, -0.00122_wp, -0.000863_wp, -0.030947_wp, &
                      -0.030562_wp, -0.019545_wp, -0.018614_wp, -0.008942_wp, -0.010571_wp, &
                      -0.023127_wp, -0.024202_wp, -0.019793_wp, -0.019387_wp, 0.036856_wp, &
                      -0.009305_wp, -0.007697_wp, 0.122222_wp, -0.016317_wp, -0.01584_wp, &
                      -0.080523_wp, -0.032986_wp, -0.032534_wp, 0.024493_wp, 0.00069_wp, &
                      0.000516_wp, 0.24774_wp, -0.008848_wp, -0.008722_wp, -0.386458_wp, &
                      -0.043126_wp, -0.042007_wp, 0.08765_wp, 0.081426_wp, 0.081723_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (127)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036438_wp, 0.018246_wp, 0.018252_wp, -0.034149_wp, 0.019122_wp, &
                      0.019114_wp, -0.070345_wp, 0.02145_wp, 0.021454_wp, 0.023294_wp, &
                      -0.036438_wp, 0.018246_wp, 0.018252_wp, -0.034149_wp, 0.019119_wp, &
                      0.019117_wp, -0.070345_wp, 0.021452_wp, 0.021451_wp, 0.023294_wp]
         case ("resp")
            charges = [ &
                      -0.0628_wp, 0.008377_wp, 0.007666_wp, 0.251211_wp, -0.041746_wp, &
                      -0.041595_wp, -0.350173_wp, 0.075642_wp, 0.076489_wp, 0.077918_wp, &
                      -0.064494_wp, 0.006691_wp, 0.007544_wp, 0.255416_wp, -0.043026_wp, &
                      -0.044068_wp, -0.344223_wp, 0.074647_wp, 0.074755_wp, 0.075771_wp]
         case ("mbis")
            charges = [ &
                      -0.215397_wp, 0.105928_wp, 0.105893_wp, -0.158845_wp, 0.099701_wp, &
                      0.099708_wp, -0.421732_wp, 0.127522_wp, 0.12755_wp, 0.129672_wp, &
                      -0.215398_wp, 0.105889_wp, 0.105933_wp, -0.158846_wp, 0.099713_wp, &
                      0.099696_wp, -0.421732_wp, 0.127542_wp, 0.12753_wp, 0.129671_wp]
         case ("chelpg")
            charges = [ &
                      -0.070173_wp, 0.009685_wp, 0.008947_wp, 0.268757_wp, -0.044781_wp, &
                      -0.044657_wp, -0.36998_wp, 0.080086_wp, 0.080917_wp, 0.082189_wp, &
                      -0.071905_wp, 0.007984_wp, 0.008817_wp, 0.273171_wp, -0.046129_wp, &
                      -0.047211_wp, -0.363868_wp, 0.079039_wp, 0.079132_wp, 0.07998_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036844_wp, 0.018599_wp, 0.018605_wp, -0.034585_wp, 0.019457_wp, &
                      0.01945_wp, -0.071608_wp, 0.021704_wp, 0.021706_wp, 0.023514_wp, &
                      -0.036844_wp, 0.018601_wp, 0.018603_wp, -0.034585_wp, 0.019455_wp, &
                      0.019453_wp, -0.071608_wp, 0.021705_wp, 0.021705_wp, 0.023514_wp]
         case ("resp")
            charges = [ &
                      -0.062689_wp, 0.007736_wp, 0.007043_wp, 0.26177_wp, -0.044033_wp, &
                      -0.043889_wp, -0.360209_wp, 0.077409_wp, 0.078249_wp, 0.079668_wp, &
                      -0.064804_wp, 0.006159_wp, 0.00702_wp, 0.265992_wp, -0.045261_wp, &
                      -0.046321_wp, -0.354385_wp, 0.07645_wp, 0.076551_wp, 0.077545_wp]
         case ("mbis")
            charges = [ &
                      -0.217336_wp, 0.106988_wp, 0.106953_wp, -0.159099_wp, 0.100551_wp, &
                      0.100558_wp, -0.428564_wp, 0.129276_wp, 0.129299_wp, 0.131373_wp, &
                      -0.217336_wp, 0.106952_wp, 0.10699_wp, -0.159099_wp, 0.100563_wp, &
                      0.100547_wp, -0.428564_wp, 0.129292_wp, 0.129282_wp, 0.131373_wp]
         case ("chelpg")
            charges = [ &
                      -0.070069_wp, 0.009037_wp, 0.008317_wp, 0.27945_wp, -0.047097_wp, &
                      -0.046981_wp, -0.380147_wp, 0.081881_wp, 0.082704_wp, 0.083967_wp, &
                      -0.072262_wp, 0.007454_wp, 0.008295_wp, 0.283894_wp, -0.048395_wp, &
                      -0.049494_wp, -0.374162_wp, 0.080871_wp, 0.080957_wp, 0.081781_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037298_wp, 0.019042_wp, 0.019048_wp, -0.035064_wp, 0.019871_wp, &
                      0.019865_wp, -0.073125_wp, 0.021977_wp, 0.021976_wp, 0.023707_wp, &
                      -0.037298_wp, 0.019047_wp, 0.019043_wp, -0.035064_wp, 0.019868_wp, &
                      0.019868_wp, -0.073125_wp, 0.021977_wp, 0.021976_wp, 0.023707_wp]
         case ("resp")
            charges = [ &
                      -0.062127_wp, 0.006776_wp, 0.006106_wp, 0.27524_wp, -0.047014_wp, &
                      -0.046878_wp, -0.372228_wp, 0.079416_wp, 0.080246_wp, 0.081601_wp, &
                      -0.064771_wp, 0.005337_wp, 0.006207_wp, 0.279486_wp, -0.048179_wp, &
                      -0.049258_wp, -0.366561_wp, 0.078501_wp, 0.078595_wp, 0.079507_wp]
         case ("mbis")
            charges = [ &
                      -0.219583_wp, 0.108262_wp, 0.108226_wp, -0.159299_wp, 0.101562_wp, &
                      0.101571_wp, -0.436664_wp, 0.131315_wp, 0.131331_wp, 0.133279_wp, &
                      -0.219583_wp, 0.108231_wp, 0.108258_wp, -0.1593_wp, 0.101573_wp, &
                      0.10156_wp, -0.436664_wp, 0.131328_wp, 0.131318_wp, 0.133279_wp]
         case ("chelpg")
            charges = [ &
                      -0.069475_wp, 0.008057_wp, 0.007359_wp, 0.293066_wp, -0.050114_wp, &
                      -0.050005_wp, -0.392312_wp, 0.083918_wp, 0.084731_wp, 0.085932_wp, &
                      -0.07225_wp, 0.00662_wp, 0.007471_wp, 0.297552_wp, -0.051351_wp, &
                      -0.052469_wp, -0.386486_wp, 0.082953_wp, 0.083031_wp, 0.083772_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (128)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07035_wp, -0.033071_wp, -0.035891_wp, -0.034092_wp, -0.035402_wp, &
                      0.23323_wp, -0.298328_wp, -0.188172_wp, 0.022073_wp, 0.024674_wp, &
                      0.022083_wp, 0.021764_wp, 0.021814_wp, 0.019384_wp, 0.019337_wp, &
                      0.026154_wp, 0.026295_wp, 0.040584_wp, 0.040131_wp, 0.177768_wp]
         case ("resp")
            charges = [ &
                      -0.356858_wp, 0.216912_wp, -0.050459_wp, -0.092913_wp, -0.052985_wp, &
                      0.676737_wp, -0.562629_wp, -0.622225_wp, 0.079139_wp, 0.084959_wp, &
                      0.079464_wp, -0.026973_wp, -0.02862_wp, 0.0166_wp, 0.016382_wp, &
                      0.040062_wp, 0.043309_wp, 0.052796_wp, 0.054648_wp, 0.432653_wp]
         case ("mbis")
            charges = [ &
                      -0.420845_wp, -0.162335_wp, -0.217497_wp, -0.175922_wp, -0.448899_wp, &
                      0.843206_wp, -0.605057_wp, -0.657263_wp, 0.127841_wp, 0.131479_wp, &
                      0.127824_wp, 0.103396_wp, 0.103514_wp, 0.107273_wp, 0.10716_wp, &
                      0.120629_wp, 0.12128_wp, 0.166314_wp, 0.16633_wp, 0.461554_wp]
         case ("chelpg")
            charges = [ &
                      -0.37405_wp, 0.23114_wp, -0.056639_wp, -0.097305_wp, -0.057326_wp, &
                      0.680916_wp, -0.563891_wp, -0.624088_wp, 0.082968_wp, 0.088907_wp, &
                      0.083287_wp, -0.029555_wp, -0.031212_wp, 0.01818_wp, 0.017946_wp, &
                      0.041775_wp, 0.045026_wp, 0.054438_wp, 0.056319_wp, 0.433165_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071949_wp, -0.033565_wp, -0.035057_wp, -0.034922_wp, -0.034888_wp, &
                      0.237375_wp, -0.321459_wp, -0.189019_wp, 0.02243_wp, 0.02451_wp, &
                      0.022439_wp, 0.021929_wp, 0.021961_wp, 0.021388_wp, 0.021358_wp, &
                      0.024658_wp, 0.024621_wp, 0.045301_wp, 0.04473_wp, 0.18814_wp]
         case ("resp")
            charges = [ &
                      -0.366843_wp, 0.216182_wp, -0.028373_wp, -0.075111_wp, -0.073324_wp, &
                      0.723029_wp, -0.611151_wp, -0.648974_wp, 0.08044_wp, 0.087012_wp, &
                      0.080777_wp, -0.028038_wp, -0.029755_wp, 0.014156_wp, 0.013947_wp, &
                      0.030074_wp, 0.033082_wp, 0.060981_wp, 0.062733_wp, 0.459155_wp]
         case ("mbis")
            charges = [ &
                      -0.429472_wp, -0.163347_wp, -0.21771_wp, -0.171879_wp, -0.467138_wp, &
                      0.873541_wp, -0.643598_wp, -0.676309_wp, 0.130149_wp, 0.132922_wp, &
                      0.130133_wp, 0.10406_wp, 0.104146_wp, 0.110601_wp, 0.11052_wp, &
                      0.117799_wp, 0.118129_wp, 0.176735_wp, 0.176546_wp, 0.484153_wp]
         case ("chelpg")
            charges = [ &
                      -0.383546_wp, 0.229696_wp, -0.033533_wp, -0.078129_wp, -0.079911_wp, &
                      0.727937_wp, -0.612506_wp, -0.65078_wp, 0.084158_wp, 0.090896_wp, &
                      0.084489_wp, -0.030569_wp, -0.032296_wp, 0.015429_wp, 0.015203_wp, &
                      0.031471_wp, 0.034484_wp, 0.063086_wp, 0.064867_wp, 0.459552_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072958_wp, -0.033968_wp, -0.034208_wp, -0.035635_wp, -0.034325_wp, &
                      0.240204_wp, -0.339046_wp, -0.188453_wp, 0.02281_wp, 0.024157_wp, &
                      0.022812_wp, 0.021758_wp, 0.021763_wp, 0.02329_wp, 0.023285_wp, &
                      0.022911_wp, 0.022718_wp, 0.049154_wp, 0.048513_wp, 0.1952_wp]
         case ("resp")
            charges = [ &
                      -0.372746_wp, 0.213658_wp, -0.008234_wp, -0.064374_wp, -0.086529_wp, &
                      0.755302_wp, -0.647506_wp, -0.665343_wp, 0.081308_wp, 0.087906_wp, &
                      0.08164_wp, -0.028994_wp, -0.030789_wp, 0.012444_wp, 0.012261_wp, &
                      0.021808_wp, 0.024583_wp, 0.067408_wp, 0.069112_wp, 0.477086_wp]
         case ("mbis")
            charges = [ &
                      -0.435013_wp, -0.163937_wp, -0.217881_wp, -0.167971_wp, -0.480384_wp, &
                      0.894604_wp, -0.672498_wp, -0.68695_wp, 0.131902_wp, 0.133402_wp, &
                      0.131875_wp, 0.104018_wp, 0.104054_wp, 0.113738_wp, 0.113702_wp, &
                      0.114468_wp, 0.114511_wp, 0.184937_wp, 0.184625_wp, 0.498781_wp]
         case ("chelpg")
            charges = [ &
                      -0.389105_wp, 0.226443_wp, -0.011988_wp, -0.066871_wp, -0.094319_wp, &
                      0.760625_wp, -0.648895_wp, -0.667116_wp, 0.084945_wp, 0.091765_wp, &
                      0.08527_wp, -0.03146_wp, -0.033263_wp, 0.013336_wp, 0.013137_wp, &
                      0.023032_wp, 0.025811_wp, 0.069754_wp, 0.071487_wp, 0.47741_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (129)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.027725_wp, 0.041179_wp, -0.032243_wp, 0.043359_wp, 0.012398_wp, &
                      -0.069137_wp, -0.032243_wp, 0.043362_wp, -0.027725_wp, 0.041179_wp, &
                      -0.032433_wp, 0.040029_wp]
         case ("resp")
            charges = [ &
                      -0.170883_wp, 0.137216_wp, -0.018702_wp, 0.106357_wp, -0.058018_wp, &
                      -0.080172_wp, -0.015815_wp, 0.105773_wp, -0.175567_wp, 0.139253_wp, &
                      -0.097691_wp, 0.128249_wp]
         case ("mbis")
            charges = [ &
                      -0.12752_wp, 0.140277_wp, -0.132526_wp, 0.15852_wp, -0.009844_wp, &
                      -0.071911_wp, -0.132524_wp, 0.158519_wp, -0.127524_wp, 0.140278_wp, &
                      -0.136065_wp, 0.140319_wp]
         case ("chelpg")
            charges = [ &
                      -0.17848_wp, 0.138789_wp, -0.011863_wp, 0.105158_wp, -0.063558_wp, &
                      -0.079805_wp, -0.008916_wp, 0.104551_wp, -0.183191_wp, 0.140828_wp, &
                      -0.090637_wp, 0.127123_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.030052_wp, 0.048525_wp, -0.036013_wp, 0.04892_wp, 0.005692_wp, &
                      -0.081312_wp, -0.036012_wp, 0.048924_wp, -0.030052_wp, 0.048525_wp, &
                      -0.034298_wp, 0.047151_wp]
         case ("resp")
            charges = [ &
                      -0.183151_wp, 0.153787_wp, -0.0254_wp, 0.119825_wp, -0.075818_wp, &
                      -0.09046_wp, -0.022575_wp, 0.11927_wp, -0.187753_wp, 0.155779_wp, &
                      -0.10746_wp, 0.143955_wp]
         case ("mbis")
            charges = [ &
                      -0.137251_wp, 0.155783_wp, -0.14417_wp, 0.171662_wp, -0.017627_wp, &
                      -0.083958_wp, -0.144168_wp, 0.171663_wp, -0.137255_wp, 0.155784_wp, &
                      -0.145836_wp, 0.15537_wp]
         case ("chelpg")
            charges = [ &
                      -0.190344_wp, 0.155282_wp, -0.018927_wp, 0.1187_wp, -0.081073_wp, &
                      -0.090112_wp, -0.016048_wp, 0.118123_wp, -0.194967_wp, 0.157275_wp, &
                      -0.100819_wp, 0.142909_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.030548_wp, 0.051057_wp, -0.037396_wp, 0.050579_wp, 0.003245_wp, &
                      -0.085803_wp, -0.037395_wp, 0.050582_wp, -0.030548_wp, 0.051057_wp, &
                      -0.03451_wp, 0.04968_wp]
         case ("resp")
            charges = [ &
                      -0.18646_wp, 0.159262_wp, -0.028379_wp, 0.124204_wp, -0.081593_wp, &
                      -0.094458_wp, -0.025576_wp, 0.123658_wp, -0.191029_wp, 0.161238_wp, &
                      -0.110076_wp, 0.149209_wp]
         case ("mbis")
            charges = [ &
                      -0.139914_wp, 0.160908_wp, -0.148212_wp, 0.175734_wp, -0.020694_wp, &
                      -0.08833_wp, -0.148209_wp, 0.175734_wp, -0.139918_wp, 0.160909_wp, &
                      -0.148414_wp, 0.160405_wp]
         case ("chelpg")
            charges = [ &
                      -0.193408_wp, 0.160712_wp, -0.022154_wp, 0.123127_wp, -0.086637_wp, &
                      -0.094129_wp, -0.019299_wp, 0.122562_wp, -0.197997_wp, 0.162688_wp, &
                      -0.103671_wp, 0.148205_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (130)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.006052_wp, -0.125978_wp, 0.041118_wp, 0.041117_wp, -0.064159_wp, &
                      0.033306_wp, 0.034273_wp, 0.034273_wp]
         case ("resp")
            charges = [ &
                      -0.099538_wp, -0.167565_wp, 0.105744_wp, 0.105371_wp, -0.143162_wp, &
                      0.062346_wp, 0.068822_wp, 0.067982_wp]
         case ("mbis")
            charges = [ &
                      -0.232968_wp, -0.141868_wp, 0.163661_wp, 0.163659_wp, -0.400848_wp, &
                      0.139598_wp, 0.154383_wp, 0.154386_wp]
         case ("chelpg")
            charges = [ &
                      -0.100309_wp, -0.16737_wp, 0.106481_wp, 0.106104_wp, -0.149682_wp, &
                      0.064098_wp, 0.070753_wp, 0.069925_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.015473_wp, -0.163061_wp, 0.049469_wp, 0.049468_wp, -0.063035_wp, &
                      0.041123_wp, 0.035283_wp, 0.035283_wp]
         case ("resp")
            charges = [ &
                      -0.050759_wp, -0.225742_wp, 0.110154_wp, 0.109755_wp, -0.164349_wp, &
                      0.077062_wp, 0.07236_wp, 0.071518_wp]
         case ("mbis")
            charges = [ &
                      -0.213209_wp, -0.189101_wp, 0.174337_wp, 0.174335_wp, -0.414886_wp, &
                      0.154162_wp, 0.157181_wp, 0.157185_wp]
         case ("chelpg")
            charges = [ &
                      -0.050275_wp, -0.225687_wp, 0.110569_wp, 0.110164_wp, -0.172144_wp, &
                      0.079093_wp, 0.074553_wp, 0.073726_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.017539_wp, -0.170713_wp, 0.051151_wp, 0.05115_wp, -0.062699_wp, &
                      0.042829_wp, 0.035373_wp, 0.035374_wp]
         case ("resp")
            charges = [ &
                      -0.038615_wp, -0.238316_wp, 0.110649_wp, 0.110243_wp, -0.169912_wp, &
                      0.08044_wp, 0.073174_wp, 0.072336_wp]
         case ("mbis")
            charges = [ &
                      -0.208756_wp, -0.198882_wp, 0.176376_wp, 0.176375_wp, -0.417237_wp, &
                      0.157234_wp, 0.157444_wp, 0.157449_wp]
         case ("chelpg")
            charges = [ &
                      -0.037697_wp, -0.23831_wp, 0.110949_wp, 0.110537_wp, -0.178101_wp, &
                      0.082555_wp, 0.075445_wp, 0.074622_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (131)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.034487_wp, 0.045901_wp, -0.103308_wp, 0.045951_wp, 0.045944_wp]
         case ("resp")
            charges = [ &
                      -0.481959_wp, 0.192494_wp, -0.09633_wp, 0.192917_wp, 0.192878_wp]
         case ("mbis")
            charges = [ &
                      -0.483053_wp, 0.191119_wp, -0.090283_wp, 0.191119_wp, 0.191099_wp]
         case ("chelpg")
            charges = [ &
                      -0.489563_wp, 0.194703_wp, -0.095346_wp, 0.195124_wp, 0.195082_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.02598_wp, 0.053493_wp, -0.134674_wp, 0.053584_wp, 0.053578_wp]
         case ("resp")
            charges = [ &
                      -0.453059_wp, 0.198813_wp, -0.14438_wp, 0.199305_wp, 0.199321_wp]
         case ("mbis")
            charges = [ &
                      -0.473949_wp, 0.201343_wp, -0.130219_wp, 0.201422_wp, 0.201402_wp]
         case ("chelpg")
            charges = [ &
                      -0.460596_wp, 0.201004_wp, -0.143408_wp, 0.201493_wp, 0.201507_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.023827_wp, 0.055168_wp, -0.141872_wp, 0.055269_wp, 0.055263_wp]
         case ("resp")
            charges = [ &
                      -0.442627_wp, 0.19926_wp, -0.156204_wp, 0.199769_wp, 0.199802_wp]
         case ("mbis")
            charges = [ &
                      -0.470882_wp, 0.203396_wp, -0.13948_wp, 0.203493_wp, 0.203474_wp]
         case ("chelpg")
            charges = [ &
                      -0.450149_wp, 0.201447_wp, -0.155235_wp, 0.201953_wp, 0.201983_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (132)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.041764_wp, -0.037387_wp, -0.042167_wp, 0.008338_wp, -0.039162_wp, &
                      -0.036848_wp, -0.000446_wp, -0.070649_wp, -0.070652_wp, 0.034298_wp, &
                      0.034787_wp, 0.030663_wp, 0.032042_wp, 0.034553_wp, 0.022973_wp, &
                      0.023788_wp, 0.024698_wp, 0.022234_wp, 0.023774_wp, 0.024691_wp, &
                      0.022234_wp]
         case ("resp")
            charges = [ &
                      -0.192882_wp, -0.046452_wp, -0.306383_wp, 0.170931_wp, -0.244165_wp, &
                      -0.090175_wp, 0.384302_wp, -0.49831_wp, -0.491221_wp, 0.133097_wp, &
                      0.115151_wp, 0.135501_wp, 0.138601_wp, 0.12147_wp, -0.019161_wp, &
                      0.123388_wp, 0.105297_wp, 0.118973_wp, 0.121734_wp, 0.102958_wp, &
                      0.117345_wp]
         case ("mbis")
            charges = [ &
                      -0.15321_wp, -0.118816_wp, -0.200954_wp, 0.084197_wp, -0.173668_wp, &
                      -0.130912_wp, 0.002189_wp, -0.448316_wp, -0.448322_wp, 0.135727_wp, &
                      0.132134_wp, 0.131664_wp, 0.133156_wp, 0.133284_wp, 0.102933_wp, &
                      0.134052_wp, 0.140736_wp, 0.134636_wp, 0.134037_wp, 0.140727_wp, &
                      0.134723_wp]
         case ("chelpg")
            charges = [ &
                      -0.211174_wp, -0.026045_wp, -0.330371_wp, 0.191049_wp, -0.261471_wp, &
                      -0.073391_wp, 0.396075_wp, -0.518072_wp, -0.510963_wp, 0.136418_wp, &
                      0.111691_wp, 0.140028_wp, 0.140727_wp, 0.118836_wp, -0.020815_wp, &
                      0.128303_wp, 0.109415_wp, 0.12384_wp, 0.126629_wp, 0.10709_wp, &
                      0.1222_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.046922_wp, -0.042207_wp, -0.046076_wp, 0.005812_wp, -0.042931_wp, &
                      -0.041651_wp, -0.000317_wp, -0.071826_wp, -0.071827_wp, 0.037824_wp, &
                      0.038511_wp, 0.034856_wp, 0.036626_wp, 0.038298_wp, 0.025801_wp, &
                      0.02549_wp, 0.02402_wp, 0.023518_wp, 0.02548_wp, 0.024018_wp, &
                      0.0235_wp]
         case ("resp")
            charges = [ &
                      -0.208507_wp, -0.055046_wp, -0.321932_wp, 0.168641_wp, -0.26286_wp, &
                      -0.097932_wp, 0.404916_wp, -0.517272_wp, -0.510046_wp, 0.144115_wp, &
                      0.125279_wp, 0.147128_wp, 0.152783_wp, 0.131486_wp, -0.01666_wp, &
                      0.129945_wp, 0.106272_wp, 0.1246_wp, 0.128276_wp, 0.103913_wp, &
                      0.122902_wp]
         case ("mbis")
            charges = [ &
                      -0.165397_wp, -0.129707_wp, -0.213636_wp, 0.08424_wp, -0.186445_wp, &
                      -0.141838_wp, 0.005228_wp, -0.458057_wp, -0.458038_wp, 0.145668_wp, &
                      0.142177_wp, 0.14253_wp, 0.14494_wp, 0.143354_wp, 0.108077_wp, &
                      0.138729_wp, 0.141083_wp, 0.138622_wp, 0.138717_wp, 0.14108_wp, &
                      0.138672_wp]
         case ("chelpg")
            charges = [ &
                      -0.226253_wp, -0.035242_wp, -0.34539_wp, 0.188358_wp, -0.279727_wp, &
                      -0.081639_wp, 0.416706_wp, -0.537135_wp, -0.529894_wp, 0.147334_wp, &
                      0.121947_wp, 0.151557_wp, 0.154839_wp, 0.128942_wp, -0.018285_wp, &
                      0.134892_wp, 0.110419_wp, 0.129499_wp, 0.133204_wp, 0.108076_wp, &
                      0.12779_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.051443_wp, -0.046365_wp, -0.049336_wp, 0.003543_wp, -0.046089_wp, &
                      -0.045771_wp, -0.000166_wp, -0.07265_wp, -0.072649_wp, 0.040641_wp, &
                      0.041561_wp, 0.038519_wp, 0.040611_wp, 0.041391_wp, 0.028356_wp, &
                      0.02724_wp, 0.023054_wp, 0.024648_wp, 0.027234_wp, 0.023057_wp, &
                      0.024611_wp]
         case ("resp")
            charges = [ &
                      -0.222155_wp, -0.062428_wp, -0.334178_wp, 0.164803_wp, -0.278839_wp, &
                      -0.104029_wp, 0.422606_wp, -0.531731_wp, -0.524396_wp, 0.153269_wp, &
                      0.133697_wp, 0.157036_wp, 0.165072_wp, 0.139768_wp, -0.014372_wp, &
                      0.13573_wp, 0.105977_wp, 0.129136_wp, 0.134055_wp, 0.103605_wp, &
                      0.127372_wp]
         case ("mbis")
            charges = [ &
                      -0.17594_wp, -0.138933_wp, -0.224115_wp, 0.083727_wp, -0.197132_wp, &
                      -0.150997_wp, 0.007455_wp, -0.465349_wp, -0.465309_wp, 0.153848_wp, &
                      0.150497_wp, 0.151858_wp, 0.155046_wp, 0.15172_wp, 0.112736_wp, &
                      0.143064_wp, 0.140467_wp, 0.141907_wp, 0.143057_wp, 0.140471_wp, &
                      0.141921_wp]
         case ("chelpg")
            charges = [ &
                      -0.239476_wp, -0.043112_wp, -0.357203_wp, 0.184201_wp, -0.295346_wp, &
                      -0.08813_wp, 0.434393_wp, -0.551664_wp, -0.544318_wp, 0.156408_wp, &
                      0.130469_wp, 0.161384_wp, 0.167072_wp, 0.137296_wp, -0.015967_wp, &
                      0.140701_wp, 0.110146_wp, 0.134061_wp, 0.139008_wp, 0.10779_wp, &
                      0.132287_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (133)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068889_wp, 0.003943_wp, -0.044297_wp, -0.042755_wp, 0.003334_wp, &
                      -0.001157_wp, -0.071497_wp, -0.07151_wp, -0.03908_wp, -0.042835_wp, &
                      0.030077_wp, 0.029874_wp, 0.029719_wp, 0.030498_wp, 0.029926_wp, &
                      0.022826_wp, 0.022056_wp, 0.0233_wp, 0.024276_wp, 0.024246_wp, &
                      0.023313_wp, 0.022046_wp, 0.031445_wp, 0.031147_wp]
         case ("resp")
            charges = [ &
                      -0.53729_wp, 0.3529_wp, -0.247644_wp, -0.263725_wp, 0.082394_wp, &
                      0.42962_wp, -0.49199_wp, -0.482675_wp, -0.197355_wp, -0.2915_wp, &
                      0.14309_wp, 0.142197_wp, 0.137263_wp, 0.151713_wp, 0.140434_wp, &
                      -0.031493_wp, 0.113695_wp, 0.116756_wp, 0.104621_wp, 0.101506_wp, &
                      0.11483_wp, 0.112195_wp, 0.144435_wp, 0.156022_wp]
         case ("mbis")
            charges = [ &
                      -0.45867_wp, 0.138063_wp, -0.186272_wp, -0.189074_wp, 0.067237_wp, &
                      0.006343_wp, -0.448258_wp, -0.448418_wp, -0.160608_wp, -0.198531_wp, &
                      0.144124_wp, 0.143747_wp, 0.143286_wp, 0.133009_wp, 0.129707_wp, &
                      0.101253_wp, 0.13423_wp, 0.132964_wp, 0.140178_wp, 0.140183_wp, &
                      0.133058_wp, 0.134349_wp, 0.131019_wp, 0.137088_wp]
         case ("chelpg")
            charges = [ &
                      -0.562623_wp, 0.370879_wp, -0.252972_wp, -0.265973_wp, 0.079127_wp, &
                      0.446865_wp, -0.512799_wp, -0.503451_wp, -0.191123_wp, -0.303071_wp, &
                      0.148822_wp, 0.148012_wp, 0.142936_wp, 0.152905_wp, 0.141711_wp, &
                      -0.03374_wp, 0.118708_wp, 0.121683_wp, 0.109249_wp, 0.106139_wp, &
                      0.119722_wp, 0.117219_wp, 0.142888_wp, 0.158887_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069896_wp, 0.000614_wp, -0.048854_wp, -0.047016_wp, 0.000353_wp, &
                      -0.001132_wp, -0.072677_wp, -0.072691_wp, -0.043031_wp, -0.047028_wp, &
                      0.031695_wp, 0.031817_wp, 0.031814_wp, 0.033987_wp, 0.033439_wp, &
                      0.025284_wp, 0.023262_wp, 0.024806_wp, 0.023544_wp, 0.023513_wp, &
                      0.024819_wp, 0.023254_wp, 0.035474_wp, 0.034658_wp]
         case ("resp")
            charges = [ &
                      -0.555154_wp, 0.35206_wp, -0.258248_wp, -0.281106_wp, 0.08064_wp, &
                      0.450271_wp, -0.510981_wp, -0.500699_wp, -0.217861_wp, -0.300009_wp, &
                      0.150372_wp, 0.149238_wp, 0.144379_wp, 0.16232_wp, 0.151446_wp, &
                      -0.029701_wp, 0.119323_wp, 0.123058_wp, 0.105597_wp, 0.102355_wp, &
                      0.120769_wp, 0.117662_wp, 0.157973_wp, 0.166295_wp]
         case ("mbis")
            charges = [ &
                      -0.468939_wp, 0.138875_wp, -0.199577_wp, -0.201552_wp, 0.066226_wp, &
                      0.009653_wp, -0.457586_wp, -0.457756_wp, -0.173167_wp, -0.210702_wp, &
                      0.148558_wp, 0.14881_wp, 0.148668_wp, 0.142935_wp, 0.139508_wp, &
                      0.105759_wp, 0.138068_wp, 0.137254_wp, 0.140412_wp, 0.140418_wp, &
                      0.13735_wp, 0.138195_wp, 0.141875_wp, 0.146721_wp]
         case ("chelpg")
            charges = [ &
                      -0.580616_wp, 0.37004_wp, -0.263524_wp, -0.283452_wp, 0.077488_wp, &
                      0.467527_wp, -0.531892_wp, -0.521575_wp, -0.21174_wp, -0.311515_wp, &
                      0.156135_wp, 0.155085_wp, 0.150084_wp, 0.163506_wp, 0.152741_wp, &
                      -0.03194_wp, 0.124364_wp, 0.128014_wp, 0.11025_wp, 0.107012_wp, &
                      0.125689_wp, 0.122713_wp, 0.156451_wp, 0.169152_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070727_wp, -0.002704_wp, -0.053234_wp, -0.051129_wp, -0.002679_wp, &
                      -0.001111_wp, -0.073623_wp, -0.073635_wp, -0.046779_wp, -0.050989_wp, &
                      0.033235_wp, 0.033689_wp, 0.033826_wp, 0.037229_wp, 0.036756_wp, &
                      0.027699_wp, 0.02446_wp, 0.026552_wp, 0.022425_wp, 0.022403_wp, &
                      0.026563_wp, 0.024456_wp, 0.039337_wp, 0.037988_wp]
         case ("resp")
            charges = [ &
                      -0.56947_wp, 0.347582_wp, -0.266869_wp, -0.297427_wp, 0.077883_wp, &
                      0.469192_wp, -0.52718_wp, -0.515692_wp, -0.238372_wp, -0.305715_wp, &
                      0.156959_wp, 0.15552_wp, 0.150773_wp, 0.171968_wp, 0.161806_wp, &
                      -0.027782_wp, 0.124458_wp, 0.129318_wp, 0.105301_wp, 0.10192_wp, &
                      0.126569_wp, 0.122595_wp, 0.171065_wp, 0.1756_wp]
         case ("mbis")
            charges = [ &
                      -0.477774_wp, 0.13879_wp, -0.211969_wp, -0.213343_wp, 0.064614_wp, &
                      0.012391_wp, -0.465382_wp, -0.465549_wp, -0.184999_wp, -0.221869_wp, &
                      0.152607_wp, 0.153502_wp, 0.153656_wp, 0.152182_wp, 0.148769_wp, &
                      0.110222_wp, 0.141627_wp, 0.141704_wp, 0.139648_wp, 0.13967_wp, &
                      0.141792_wp, 0.14176_wp, 0.152194_wp, 0.155763_wp]
         case ("chelpg")
            charges = [ &
                      -0.595036_wp, 0.365594_wp, -0.272142_wp, -0.29979_wp, 0.074742_wp, &
                      0.486488_wp, -0.548182_wp, -0.536655_wp, -0.232281_wp, -0.317205_wp, &
                      0.162747_wp, 0.161392_wp, 0.156505_wp, 0.173155_wp, 0.163105_wp, &
                      -0.030019_wp, 0.129523_wp, 0.134297_wp, 0.109976_wp, 0.106599_wp, &
                      0.131511_wp, 0.12767_wp, 0.169551_wp, 0.178455_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (134)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.055582_wp, 0.007589_wp, -0.05544_wp, 0.007558_wp, -0.055657_wp, &
                      0.007499_wp, -0.06815_wp, -0.068192_wp, -0.068152_wp, 0.029623_wp, &
                      0.027394_wp, 0.026931_wp, 0.026624_wp, 0.030012_wp, 0.029153_wp, &
                      0.030003_wp, 0.0294_wp, 0.030168_wp, 0.029893_wp, 0.029227_wp, &
                      0.030085_wp]
         case ("resp")
            charges = [ &
                      -0.624023_wp, 0.512819_wp, -0.619571_wp, 0.50911_wp, -0.623322_wp, &
                      0.507805_wp, -0.454032_wp, -0.466753_wp, -0.465069_wp, 0.120803_wp, &
                      0.211707_wp, 0.212387_wp, 0.212787_wp, 0.121946_wp, 0.119734_wp, &
                      0.123475_wp, 0.119427_wp, 0.124612_wp, 0.118844_wp, 0.115881_wp, &
                      0.121432_wp]
         case ("mbis")
            charges = [ &
                      -0.292174_wp, 0.178696_wp, -0.292451_wp, 0.178292_wp, -0.2911_wp, &
                      0.178988_wp, -0.462879_wp, -0.463068_wp, -0.462333_wp, 0.14418_wp, &
                      0.144444_wp, 0.142752_wp, 0.14136_wp, 0.144999_wp, 0.143429_wp, &
                      0.144919_wp, 0.14372_wp, 0.144953_wp, 0.144787_wp, 0.14353_wp, &
                      0.144941_wp]
         case ("chelpg")
            charges = [ &
                      -0.674092_wp, 0.565951_wp, -0.670655_wp, 0.560776_wp, -0.672117_wp, &
                      0.559899_wp, -0.480808_wp, -0.493845_wp, -0.489968_wp, 0.125402_wp, &
                      0.221242_wp, 0.221432_wp, 0.221361_wp, 0.126958_wp, 0.125024_wp, &
                      0.128602_wp, 0.124128_wp, 0.129237_wp, 0.123866_wp, 0.12113_wp, &
                      0.126477_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06039_wp, 0.003987_wp, -0.060039_wp, 0.003924_wp, -0.060546_wp, &
                      0.003833_wp, -0.069275_wp, -0.069296_wp, -0.06926_wp, 0.031701_wp, &
                      0.03094_wp, 0.030472_wp, 0.030194_wp, 0.031908_wp, 0.031366_wp, &
                      0.031851_wp, 0.031573_wp, 0.031888_wp, 0.031835_wp, 0.031439_wp, &
                      0.031879_wp]
         case ("resp")
            charges = [ &
                      -0.637719_wp, 0.509452_wp, -0.632695_wp, 0.505839_wp, -0.637374_wp, &
                      0.504557_wp, -0.469657_wp, -0.482251_wp, -0.480836_wp, 0.127614_wp, &
                      0.223507_wp, 0.224241_wp, 0.224755_wp, 0.128662_wp, 0.126496_wp, &
                      0.13056_wp, 0.126197_wp, 0.131829_wp, 0.125598_wp, 0.122673_wp, &
                      0.128552_wp]
         case ("mbis")
            charges = [ &
                      -0.309604_wp, 0.18056_wp, -0.309107_wp, 0.180189_wp, -0.30891_wp, &
                      0.181067_wp, -0.473774_wp, -0.473871_wp, -0.47317_wp, 0.149586_wp, &
                      0.154775_wp, 0.153274_wp, 0.152047_wp, 0.150033_wp, 0.149111_wp, &
                      0.149855_wp, 0.149312_wp, 0.149654_wp, 0.149935_wp, 0.149226_wp, &
                      0.149797_wp]
         case ("chelpg")
            charges = [ &
                      -0.687847_wp, 0.562666_wp, -0.683849_wp, 0.557584_wp, -0.686219_wp, &
                      0.556737_wp, -0.496585_wp, -0.509487_wp, -0.505871_wp, 0.132247_wp, &
                      0.233061_wp, 0.2333_wp, 0.23334_wp, 0.13371_wp, 0.131822_wp, &
                      0.135722_wp, 0.130933_wp, 0.136487_wp, 0.130657_wp, 0.12796_wp, &
                      0.133633_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.064898_wp, 0.000489_wp, -0.064345_wp, 0.000401_wp, -0.06512_wp, &
                      0.000277_wp, -0.070184_wp, -0.070181_wp, -0.070149_wp, 0.033654_wp, &
                      0.034224_wp, 0.033728_wp, 0.033475_wp, 0.033681_wp, 0.033468_wp, &
                      0.03359_wp, 0.033636_wp, 0.033496_wp, 0.033647_wp, 0.033538_wp, &
                      0.033556_wp]
         case ("resp")
            charges = [ &
                      -0.646836_wp, 0.501962_wp, -0.641307_wp, 0.49838_wp, -0.646786_wp, &
                      0.497128_wp, -0.481926_wp, -0.494424_wp, -0.493186_wp, 0.13364_wp, &
                      0.233771_wp, 0.234492_wp, 0.235109_wp, 0.134606_wp, 0.132543_wp, &
                      0.136971_wp, 0.132224_wp, 0.138347_wp, 0.13157_wp, 0.128737_wp, &
                      0.134986_wp]
         case ("mbis")
            charges = [ &
                      -0.325454_wp, 0.181522_wp, -0.324213_wp, 0.18118_wp, -0.325094_wp, &
                      0.18222_wp, -0.482964_wp, -0.482988_wp, -0.482307_wp, 0.154477_wp, &
                      0.164269_wp, 0.162915_wp, 0.161833_wp, 0.154558_wp, 0.154304_wp, &
                      0.154308_wp, 0.154415_wp, 0.153861_wp, 0.154561_wp, 0.154429_wp, &
                      0.154155_wp]
         case ("chelpg")
            charges = [ &
                      -0.697009_wp, 0.555252_wp, -0.69251_wp, 0.550192_wp, -0.695668_wp, &
                      0.549372_wp, -0.508969_wp, -0.521771_wp, -0.518324_wp, 0.138296_wp, &
                      0.243335_wp, 0.243559_wp, 0.2437_wp, 0.139679_wp, 0.137894_wp, &
                      0.142158_wp, 0.136984_wp, 0.143028_wp, 0.136656_wp, 0.134052_wp, &
                      0.140094_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (135)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.014882_wp, -0.251085_wp, 0.023413_wp, 0.023432_wp, 0.037317_wp, &
                      0.15204_wp]
         case ("resp")
            charges = [ &
                      0.077249_wp, -0.578264_wp, 0.018375_wp, 0.019096_wp, 0.082958_wp, &
                      0.380587_wp]
         case ("mbis")
            charges = [ &
                      -0.070712_wp, -0.597779_wp, 0.071842_wp, 0.071843_wp, 0.108992_wp, &
                      0.415813_wp]
         case ("chelpg")
            charges = [ &
                      0.080101_wp, -0.57913_wp, 0.017593_wp, 0.018315_wp, 0.082339_wp, &
                      0.380781_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.014002_wp, -0.280082_wp, 0.030533_wp, 0.030562_wp, 0.036446_wp, &
                      0.168537_wp]
         case ("resp")
            charges = [ &
                      0.108669_wp, -0.673588_wp, 0.027487_wp, 0.0282_wp, 0.073295_wp, &
                      0.435937_wp]
         case ("mbis")
            charges = [ &
                      -0.07039_wp, -0.670689_wp, 0.084785_wp, 0.084804_wp, 0.108269_wp, &
                      0.463219_wp]
         case ("chelpg")
            charges = [ &
                      0.112007_wp, -0.674552_wp, 0.026578_wp, 0.027291_wp, 0.07255_wp, &
                      0.436126_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.013983_wp, -0.281146_wp, 0.030826_wp, 0.030856_wp, 0.036355_wp, &
                      0.169124_wp]
         case ("resp")
            charges = [ &
                      0.109887_wp, -0.677152_wp, 0.027887_wp, 0.0286_wp, 0.072809_wp, &
                      0.437969_wp]
         case ("mbis")
            charges = [ &
                      -0.070319_wp, -0.67335_wp, 0.08531_wp, 0.08533_wp, 0.108124_wp, &
                      0.464903_wp]
         case ("chelpg")
            charges = [ &
                      0.113234_wp, -0.678119_wp, 0.026977_wp, 0.027689_wp, 0.07206_wp, &
                      0.438159_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (136)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.033807_wp, -0.022477_wp, -0.031998_wp, -0.017728_wp, -0.01535_wp, &
                      -0.024206_wp, 0.222309_wp, -0.291477_wp, -0.138073_wp, 0.029234_wp, &
                      0.038895_wp, 0.040566_wp, 0.039991_wp, 0.043326_wp, 0.038107_wp, &
                      0.042304_wp, 0.040192_wp, 0.040189_wp]
         case ("resp")
            charges = [ &
                      -0.084472_wp, -0.140398_wp, -0.115964_wp, -0.072372_wp, -0.126315_wp, &
                      -0.121866_wp, 0.684248_wp, -0.524047_wp, -0.26557_wp, -0.260525_wp, &
                      0.117274_wp, 0.133171_wp, 0.124257_wp, 0.11995_wp, 0.112546_wp, &
                      0.142389_wp, 0.138664_wp, 0.139031_wp]
         case ("mbis")
            charges = [ &
                      -0.159698_wp, -0.098232_wp, -0.159836_wp, -0.07885_wp, -0.172558_wp, &
                      -0.090545_wp, 0.767074_wp, -0.579238_wp, -0.379642_wp, -0.121571_wp, &
                      0.141351_wp, 0.135294_wp, 0.142627_wp, 0.153123_wp, 0.147193_wp, &
                      0.11816_wp, 0.11768_wp, 0.117665_wp]
         case ("chelpg")
            charges = [ &
                      -0.084744_wp, -0.139759_wp, -0.116752_wp, -0.070569_wp, -0.13151_wp, &
                      -0.120464_wp, 0.68998_wp, -0.525696_wp, -0.264886_wp, -0.271176_wp, &
                      0.117195_wp, 0.133114_wp, 0.124363_wp, 0.119829_wp, 0.112535_wp, &
                      0.14544_wp, 0.141372_wp, 0.14173_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.034412_wp, -0.022082_wp, -0.033786_wp, -0.024448_wp, -0.022406_wp, &
                      -0.026976_wp, 0.224411_wp, -0.32135_wp, -0.136862_wp, 0.030116_wp, &
                      0.047822_wp, 0.049943_wp, 0.048136_wp, 0.042528_wp, 0.042456_wp, &
                      0.051071_wp, 0.042918_wp, 0.042919_wp]
         case ("resp")
            charges = [ &
                      -0.091023_wp, -0.149556_wp, -0.126614_wp, -0.078746_wp, -0.155095_wp, &
                      -0.12938_wp, 0.736889_wp, -0.587526_wp, -0.285143_wp, -0.255647_wp, &
                      0.134951_wp, 0.152115_wp, 0.142145_wp, 0.124346_wp, 0.125034_wp, &
                      0.15813_wp, 0.14237_wp, 0.142751_wp]
         case ("mbis")
            charges = [ &
                      -0.169074_wp, -0.104831_wp, -0.169521_wp, -0.088748_wp, -0.192232_wp, &
                      -0.096131_wp, 0.79695_wp, -0.627003_wp, -0.383135_wp, -0.133292_wp, &
                      0.159002_wp, 0.153221_wp, 0.159371_wp, 0.156247_wp, 0.157233_wp, &
                      0.134593_wp, 0.123678_wp, 0.12367_wp]
         case ("chelpg")
            charges = [ &
                      -0.091638_wp, -0.148594_wp, -0.12771_wp, -0.076595_wp, -0.160737_wp, &
                      -0.127596_wp, 0.742811_wp, -0.589211_wp, -0.284488_wp, -0.26635_wp, &
                      0.134933_wp, 0.152002_wp, 0.142303_wp, 0.124168_wp, 0.124957_wp, &
                      0.161202_wp, 0.145086_wp, 0.145459_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.034131_wp, -0.021538_wp, -0.034035_wp, -0.02623_wp, -0.024168_wp, &
                      -0.027446_wp, 0.224556_wp, -0.32786_wp, -0.136351_wp, 0.030375_wp, &
                      0.049852_wp, 0.05209_wp, 0.049827_wp, 0.041916_wp, 0.043417_wp, &
                      0.053078_wp, 0.043321_wp, 0.043324_wp]
         case ("resp")
            charges = [ &
                      -0.091658_wp, -0.150392_wp, -0.128549_wp, -0.080801_wp, -0.161533_wp, &
                      -0.13093_wp, 0.747135_wp, -0.601087_wp, -0.288525_wp, -0.253988_wp, &
                      0.138651_wp, 0.156053_wp, 0.145784_wp, 0.124726_wp, 0.127719_wp, &
                      0.161609_wp, 0.142701_wp, 0.143085_wp]
         case ("mbis")
            charges = [ &
                      -0.170362_wp, -0.105386_wp, -0.171207_wp, -0.091343_wp, -0.196885_wp, &
                      -0.097112_wp, 0.802675_wp, -0.637233_wp, -0.383349_wp, -0.135255_wp, &
                      0.162732_wp, 0.157014_wp, 0.162751_wp, 0.15629_wp, 0.159344_wp, &
                      0.138232_wp, 0.124548_wp, 0.124542_wp]
         case ("chelpg")
            charges = [ &
                      -0.092291_wp, -0.149416_wp, -0.129654_wp, -0.078643_wp, -0.167204_wp, &
                      -0.129126_wp, 0.753086_wp, -0.60278_wp, -0.287878_wp, -0.264692_wp, &
                      0.138638_wp, 0.155937_wp, 0.145944_wp, 0.124548_wp, 0.12764_wp, &
                      0.164683_wp, 0.145417_wp, 0.145792_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (137)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.03146_wp, -0.036113_wp, 0.230451_wp, -0.291205_wp, -0.151027_wp, &
                      0.026645_wp, -0.068503_wp, 0.026624_wp, 0.026613_wp, 0.038366_wp, &
                      0.038377_wp, 0.041606_wp, 0.038543_wp, 0.038534_wp, 0.027475_wp, &
                      0.022524_wp, 0.022542_wp]
         case ("resp")
            charges = [ &
                      0.393366_wp, -0.510858_wp, 0.851959_wp, -0.578272_wp, -0.390086_wp, &
                      -0.145191_wp, -0.367928_wp, -0.047694_wp, -0.048366_wp, 0.124585_wp, &
                      0.124754_wp, 0.123033_wp, 0.110213_wp, 0.110242_wp, 0.073632_wp, &
                      0.088338_wp, 0.088271_wp]
         case ("mbis")
            charges = [ &
                      -0.118452_wp, -0.467216_wp, 0.807323_wp, -0.58146_wp, -0.398519_wp, &
                      -0.11756_wp, -0.426502_wp, 0.114471_wp, 0.114224_wp, 0.167051_wp, &
                      0.167078_wp, 0.117339_wp, 0.114433_wp, 0.114487_wp, 0.135728_wp, &
                      0.128786_wp, 0.128781_wp]
         case ("chelpg")
            charges = [ &
                      0.415584_wp, -0.537898_wp, 0.862976_wp, -0.580636_wp, -0.388891_wp, &
                      -0.156699_wp, -0.3869_wp, -0.051248_wp, -0.051947_wp, 0.130259_wp, &
                      0.130415_wp, 0.126341_wp, 0.112892_wp, 0.112909_wp, 0.077134_wp, &
                      0.092902_wp, 0.092806_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.032446_wp, -0.035108_wp, 0.235821_wp, -0.322369_wp, -0.152483_wp, &
                      0.026732_wp, -0.068597_wp, 0.024266_wp, 0.024351_wp, 0.045202_wp, &
                      0.045224_wp, 0.049033_wp, 0.040979_wp, 0.040952_wp, 0.027987_wp, &
                      0.025216_wp, 0.025231_wp]
         case ("resp")
            charges = [ &
                      0.435248_wp, -0.537933_wp, 0.913491_wp, -0.647979_wp, -0.420144_wp, &
                      -0.1331_wp, -0.37357_wp, -0.062913_wp, -0.063506_wp, 0.137083_wp, &
                      0.137278_wp, 0.134914_wp, 0.112527_wp, 0.112585_wp, 0.071305_wp, &
                      0.092397_wp, 0.092315_wp]
         case ("mbis")
            charges = [ &
                      -0.109959_wp, -0.491439_wp, 0.845223_wp, -0.633348_wp, -0.408346_wp, &
                      -0.128137_wp, -0.434801_wp, 0.109713_wp, 0.109638_wp, 0.181697_wp, &
                      0.181744_wp, 0.131597_wp, 0.120053_wp, 0.12008_wp, 0.137385_wp, &
                      0.13445_wp, 0.13444_wp]
         case ("chelpg")
            charges = [ &
                      0.457927_wp, -0.565265_wp, 0.924582_wp, -0.65035_wp, -0.418962_wp, &
                      -0.144542_wp, -0.392739_wp, -0.066579_wp, -0.067199_wp, 0.142805_wp, &
                      0.142988_wp, 0.138209_wp, 0.115184_wp, 0.11523_wp, 0.074834_wp, &
                      0.096993_wp, 0.096883_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.03276_wp, -0.034752_wp, 0.237058_wp, -0.330825_wp, -0.152376_wp, &
                      0.02684_wp, -0.068495_wp, 0.023331_wp, 0.023446_wp, 0.047193_wp, &
                      0.047218_wp, 0.051147_wp, 0.041401_wp, 0.041369_wp, 0.027977_wp, &
                      0.026105_wp, 0.026116_wp]
         case ("resp")
            charges = [ &
                      0.444641_wp, -0.543139_wp, 0.927867_wp, -0.666489_wp, -0.426613_wp, &
                      -0.129499_wp, -0.373541_wp, -0.067184_wp, -0.067748_wp, 0.140413_wp, &
                      0.140613_wp, 0.138223_wp, 0.112661_wp, 0.112728_wp, 0.070125_wp, &
                      0.093518_wp, 0.093426_wp]
         case ("mbis")
            charges = [ &
                      -0.107422_wp, -0.497583_wp, 0.854511_wp, -0.647272_wp, -0.410038_wp, &
                      -0.130368_wp, -0.436605_wp, 0.107874_wp, 0.107851_wp, 0.185798_wp, &
                      0.18585_wp, 0.135516_wp, 0.121045_wp, 0.121064_wp, 0.13745_wp, &
                      0.136169_wp, 0.136153_wp]
         case ("chelpg")
            charges = [ &
                      0.467428_wp, -0.570539_wp, 0.938978_wp, -0.668873_wp, -0.425434_wp, &
                      -0.140905_wp, -0.392739_wp, -0.070877_wp, -0.071468_wp, 0.146146_wp, &
                      0.146334_wp, 0.141506_wp, 0.11531_wp, 0.115365_wp, 0.073655_wp, &
                      0.098116_wp, 0.097996_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (138)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.231027_wp, -0.301348_wp, -0.149163_wp, -0.065052_wp, 0.046573_wp, &
                      0.048798_wp, 0.041252_wp, 0.027222_wp, 0.042009_wp, 0.03921_wp, &
                      0.039467_wp]
         case ("resp")
            charges = [ &
                      0.865386_wp, -0.569273_wp, -0.351893_wp, -0.731672_wp, 0.205903_wp, &
                      0.202441_wp, 0.19215_wp, -0.184194_wp, 0.130129_wp, 0.119826_wp, &
                      0.121196_wp]
         case ("mbis")
            charges = [ &
                      0.847426_wp, -0.596274_wp, -0.397029_wp, -0.66886_wp, 0.199076_wp, &
                      0.19427_wp, 0.191287_wp, -0.119744_wp, 0.11804_wp, 0.115541_wp, &
                      0.116262_wp]
         case ("chelpg")
            charges = [ &
                      0.873576_wp, -0.570864_wp, -0.351051_wp, -0.750481_wp, 0.210583_wp, &
                      0.206841_wp, 0.196768_wp, -0.194366_wp, 0.133091_wp, 0.12228_wp, &
                      0.123623_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.236113_wp, -0.336603_wp, -0.149878_wp, -0.064347_wp, 0.049605_wp, &
                      0.056844_wp, 0.047421_wp, 0.027333_wp, 0.049933_wp, 0.041793_wp, &
                      0.041781_wp]
         case ("resp")
            charges = [ &
                      0.914651_wp, -0.639209_wp, -0.377743_wp, -0.755441_wp, 0.214759_wp, &
                      0.218567_wp, 0.206934_wp, -0.169244_wp, 0.141443_wp, 0.12223_wp, &
                      0.123053_wp]
         case ("mbis")
            charges = [ &
                      0.886086_wp, -0.653125_wp, -0.405545_wp, -0.696176_wp, 0.207333_wp, &
                      0.211037_wp, 0.204893_wp, -0.130888_wp, 0.133172_wp, 0.121485_wp, &
                      0.121725_wp]
         case ("chelpg")
            charges = [ &
                      0.922852_wp, -0.640794_wp, -0.376923_wp, -0.774308_wp, 0.219454_wp, &
                      0.22298_wp, 0.211569_wp, -0.179337_wp, 0.144385_wp, 0.124662_wp, &
                      0.125458_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.236992_wp, -0.343922_wp, -0.149663_wp, -0.06412_wp, 0.05006_wp, &
                      0.058503_wp, 0.048825_wp, 0.027421_wp, 0.051682_wp, 0.042151_wp, &
                      0.042065_wp]
         case ("resp")
            charges = [ &
                      0.92336_wp, -0.653541_wp, -0.381933_wp, -0.759103_wp, 0.216048_wp, &
                      0.22174_wp, 0.210003_wp, -0.165825_wp, 0.14391_wp, 0.122352_wp, &
                      0.122988_wp]
         case ("mbis")
            charges = [ &
                      0.893329_wp, -0.664808_wp, -0.406612_wp, -0.70119_wp, 0.208642_wp, &
                      0.21439_wp, 0.207817_wp, -0.13272_wp, 0.136402_wp, 0.122321_wp, &
                      0.122427_wp]
         case ("chelpg")
            charges = [ &
                      0.93157_wp, -0.655137_wp, -0.381119_wp, -0.777986_wp, 0.220749_wp, &
                      0.226159_wp, 0.214642_wp, -0.175891_wp, 0.146845_wp, 0.124779_wp, &
                      0.125389_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (139)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.209528_wp, -0.286035_wp, -0.142833_wp, 0.057427_wp, 0.03141_wp, &
                      0.045469_wp, 0.042515_wp, 0.042518_wp]
         case ("resp")
            charges = [ &
                      0.600545_wp, -0.505236_wp, -0.338563_wp, 0.010817_wp, -0.112305_wp, &
                      0.111598_wp, 0.116587_wp, 0.116557_wp]
         case ("mbis")
            charges = [ &
                      0.612719_wp, -0.53776_wp, -0.371969_wp, 0.0502_wp, -0.111431_wp, &
                      0.12091_wp, 0.118663_wp, 0.118666_wp]
         case ("chelpg")
            charges = [ &
                      0.603259_wp, -0.506117_wp, -0.337776_wp, 0.009825_wp, -0.120814_wp, &
                      0.114039_wp, 0.118805_wp, 0.118781_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.217461_wp, -0.323047_wp, -0.14327_wp, 0.071757_wp, 0.03211_wp, &
                      0.053525_wp, 0.045729_wp, 0.045733_wp]
         case ("resp")
            charges = [ &
                      0.651762_wp, -0.577991_wp, -0.362478_wp, 0.026096_wp, -0.104958_wp, &
                      0.124831_wp, 0.12137_wp, 0.121369_wp]
         case ("mbis")
            charges = [ &
                      0.639288_wp, -0.595459_wp, -0.379987_wp, 0.071787_wp, -0.123035_wp, &
                      0.136199_wp, 0.125599_wp, 0.125605_wp]
         case ("chelpg")
            charges = [ &
                      0.654414_wp, -0.578851_wp, -0.361695_wp, 0.025128_wp, -0.113376_wp, &
                      0.127243_wp, 0.123567_wp, 0.123571_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.218621_wp, -0.328763_wp, -0.143024_wp, 0.073907_wp, 0.032258_wp, &
                      0.054831_wp, 0.046082_wp, 0.046087_wp]
         case ("resp")
            charges = [ &
                      0.659327_wp, -0.589241_wp, -0.365397_wp, 0.028419_wp, -0.103646_wp, &
                      0.126988_wp, 0.121773_wp, 0.121776_wp]
         case ("mbis")
            charges = [ &
                      0.643187_wp, -0.604358_wp, -0.380682_wp, 0.075018_wp, -0.124513_wp, &
                      0.138612_wp, 0.126364_wp, 0.12637_wp]
         case ("chelpg")
            charges = [ &
                      0.661972_wp, -0.590103_wp, -0.364616_wp, 0.027456_wp, -0.112042_wp, &
                      0.129395_wp, 0.123965_wp, 0.123973_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (140)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.230092_wp, -0.291374_wp, -0.150997_wp, -0.033371_wp, 0.026759_wp, &
                      -0.066207_wp, 0.039992_wp, 0.039554_wp, 0.041647_wp, 0.038619_wp, &
                      0.038632_wp, 0.029811_wp, 0.029935_wp, 0.026905_wp]
         case ("resp")
            charges = [ &
                      0.782289_wp, -0.550875_wp, -0.42573_wp, -0.139108_wp, -0.096424_wp, &
                      -0.207457_wp, 0.062079_wp, 0.061859_wp, 0.116269_wp, 0.098728_wp, &
                      0.098235_wp, 0.063073_wp, 0.062652_wp, 0.07441_wp]
         case ("mbis")
            charges = [ &
                      0.800418_wp, -0.57614_wp, -0.400074_wp, -0.410855_wp, -0.116859_wp, &
                      -0.37449_wp, 0.160791_wp, 0.160786_wp, 0.117328_wp, 0.114383_wp, &
                      0.114387_wp, 0.140591_wp, 0.141388_wp, 0.128343_wp]
         case ("chelpg")
            charges = [ &
                      0.788472_wp, -0.552125_wp, -0.426394_wp, -0.146309_wp, -0.103447_wp, &
                      -0.213601_wp, 0.064233_wp, 0.064056_wp, 0.118534_wp, 0.100435_wp, &
                      0.09995_wp, 0.064972_wp, 0.064571_wp, 0.076654_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.235567_wp, -0.323377_wp, -0.152346_wp, -0.031838_wp, 0.026838_wp, &
                      -0.067909_wp, 0.047291_wp, 0.046654_wp, 0.049229_wp, 0.041063_wp, &
                      0.041102_wp, 0.027671_wp, 0.027463_wp, 0.032587_wp]
         case ("resp")
            charges = [ &
                      0.836576_wp, -0.620283_wp, -0.457392_wp, -0.137384_wp, -0.081507_wp, &
                      -0.202885_wp, 0.070289_wp, 0.069701_wp, 0.127918_wp, 0.100609_wp, &
                      0.100218_wp, 0.054912_wp, 0.053913_wp, 0.085315_wp]
         case ("mbis")
            charges = [ &
                      0.837653_wp, -0.628848_wp, -0.409762_wp, -0.430553_wp, -0.127718_wp, &
                      -0.378489_wp, 0.175407_wp, 0.175055_wp, 0.131894_wp, 0.12006_wp, &
                      0.120108_wp, 0.137729_wp, 0.137934_wp, 0.139526_wp]
         case ("chelpg")
            charges = [ &
                      0.84283_wp, -0.621542_wp, -0.458169_wp, -0.144784_wp, -0.088119_wp, &
                      -0.208869_wp, 0.072492_wp, 0.071948_wp, 0.130081_wp, 0.102205_wp, &
                      0.101822_wp, 0.056777_wp, 0.055797_wp, 0.087529_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.236691_wp, -0.331189_wp, -0.152216_wp, -0.031378_wp, 0.026932_wp, &
                      -0.06827_wp, 0.049175_wp, 0.048501_wp, 0.051165_wp, 0.041438_wp, &
                      0.041481_wp, 0.026885_wp, 0.026585_wp, 0.034197_wp]
         case ("resp")
            charges = [ &
                      0.847858_wp, -0.636958_wp, -0.463467_wp, -0.135941_wp, -0.077724_wp, &
                      -0.201801_wp, 0.072341_wp, 0.071677_wp, 0.130849_wp, 0.100658_wp, &
                      0.100285_wp, 0.052452_wp, 0.05129_wp, 0.08848_wp]
         case ("mbis")
            charges = [ &
                      0.84581_wp, -0.641574_wp, -0.411225_wp, -0.435035_wp, -0.129778_wp, &
                      -0.378898_wp, 0.17906_wp, 0.178642_wp, 0.135487_wp, 0.120951_wp, &
                      0.121007_wp, 0.136464_wp, 0.136504_wp, 0.142581_wp]
         case ("chelpg")
            charges = [ &
                      0.854131_wp, -0.638229_wp, -0.464272_wp, -0.143385_wp, -0.08421_wp, &
                      -0.207751_wp, 0.074555_wp, 0.073936_wp, 0.132979_wp, 0.102224_wp, &
                      0.101858_wp, 0.054311_wp, 0.053169_wp, 0.090685_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (141)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.139903_wp, 0.056642_wp, -0.009657_wp, -0.06713_wp, -0.06398_wp, &
                      -0.035834_wp, -0.035623_wp, -0.062766_wp, 0.113533_wp, 0.037649_wp, &
                      0.022566_wp, 0.03195_wp, 0.026784_wp, 0.028282_wp, 0.034266_wp, &
                      0.033907_wp, 0.029315_wp]
         case ("resp")
            charges = [ &
                      -0.564119_wp, 0.465302_wp, -0.119437_wp, -0.239083_wp, -0.429829_wp, &
                      -0.007986_wp, -0.13785_wp, -0.235476_wp, 0.332021_wp, 0.099599_wp, &
                      0.059284_wp, 0.088097_wp, 0.134903_wp, 0.162851_wp, 0.116196_wp, &
                      0.138762_wp, 0.136766_wp]
         case ("mbis")
            charges = [ &
                      -0.534433_wp, 0.363964_wp, -0.168549_wp, -0.322712_wp, -0.320228_wp, &
                      -0.06411_wp, -0.067539_wp, -0.228811_wp, 0.33665_wp, 0.115227_wp, &
                      0.089465_wp, 0.108907_wp, 0.153812_wp, 0.143784_wp, 0.126558_wp, &
                      0.12689_wp, 0.141126_wp]
         case ("chelpg")
            charges = [ &
                      -0.570941_wp, 0.494887_wp, -0.126817_wp, -0.263565_wp, -0.457006_wp, &
                      0.015038_wp, -0.117966_wp, -0.25647_wp, 0.333265_wp, 0.102346_wp, &
                      0.060826_wp, 0.090042_wp, 0.139787_wp, 0.16794_wp, 0.11227_wp, &
                      0.135742_wp, 0.14062_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.142052_wp, 0.050424_wp, -0.009591_wp, -0.075091_wp, -0.07022_wp, &
                      -0.044983_wp, -0.044992_wp, -0.074694_wp, 0.12678_wp, 0.04437_wp, &
                      0.025316_wp, 0.033743_wp, 0.032857_wp, 0.037302_wp, 0.039422_wp, &
                      0.038909_wp, 0.032501_wp]
         case ("resp")
            charges = [ &
                      -0.648758_wp, 0.5154_wp, -0.081817_wp, -0.28166_wp, -0.470496_wp, &
                      -0.019574_wp, -0.151174_wp, -0.273095_wp, 0.380698_wp, 0.103888_wp, &
                      0.059985_wp, 0.080767_wp, 0.157533_wp, 0.186541_wp, 0.132424_wp, &
                      0.155359_wp, 0.153978_wp]
         case ("mbis")
            charges = [ &
                      -0.559178_wp, 0.36808_wp, -0.175284_wp, -0.346623_wp, -0.342715_wp, &
                      -0.081463_wp, -0.08435_wp, -0.255633_wp, 0.366_wp, 0.127859_wp, &
                      0.094911_wp, 0.112739_wp, 0.171319_wp, 0.165815_wp, 0.141781_wp, &
                      0.141914_wp, 0.154831_wp]
         case ("chelpg")
            charges = [ &
                      -0.655993_wp, 0.544685_wp, -0.087857_wp, -0.305769_wp, -0.496984_wp, &
                      0.002532_wp, -0.131858_wp, -0.293319_wp, 0.381968_wp, 0.106313_wp, &
                      0.061175_wp, 0.082382_wp, 0.162374_wp, 0.191522_wp, 0.128708_wp, &
                      0.152438_wp, 0.157684_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.142256_wp, 0.049041_wp, -0.009432_wp, -0.077014_wp, -0.071485_wp, &
                      -0.047247_wp, -0.047395_wp, -0.07781_wp, 0.130043_wp, 0.046119_wp, &
                      0.026039_wp, 0.033995_wp, 0.034264_wp, 0.039616_wp, 0.040549_wp, &
                      0.039947_wp, 0.033028_wp]
         case ("resp")
            charges = [ &
                      -0.670097_wp, 0.527955_wp, -0.070356_wp, -0.292042_wp, -0.479402_wp, &
                      -0.022637_wp, -0.154458_wp, -0.282723_wp, 0.392818_wp, 0.104623_wp, &
                      0.059981_wp, 0.077953_wp, 0.162853_wp, 0.192209_wp, 0.136231_wp, &
                      0.159175_wp, 0.157915_wp]
         case ("mbis")
            charges = [ &
                      -0.564394_wp, 0.36888_wp, -0.176481_wp, -0.352267_wp, -0.34739_wp, &
                      -0.085753_wp, -0.088468_wp, -0.262461_wp, 0.372888_wp, 0.131028_wp, &
                      0.096233_wp, 0.113219_wp, 0.175408_wp, 0.171172_wp, 0.145262_wp, &
                      0.145292_wp, 0.157834_wp]
         case ("chelpg")
            charges = [ &
                      -0.677491_wp, 0.557195_wp, -0.075845_wp, -0.316078_wp, -0.505728_wp, &
                      -0.000765_wp, -0.135282_wp, -0.302763_wp, 0.3941_wp, 0.106912_wp, &
                      0.061029_wp, 0.079438_wp, 0.167681_wp, 0.197162_wp, 0.132573_wp, &
                      0.156283_wp, 0.161578_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (142)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.039665_wp, -0.006025_wp, -0.039664_wp, -0.038372_wp, -0.0386_wp, &
                      -0.038373_wp, -0.071875_wp, 0.019981_wp, 0.017123_wp, 0.016799_wp, &
                      0.017121_wp, 0.019982_wp, 0.021153_wp, 0.017707_wp, 0.017784_wp, &
                      0.021114_wp, 0.017708_wp, 0.021152_wp, 0.020663_wp, 0.022143_wp, &
                      0.022143_wp]
         case ("resp")
            charges = [ &
                      -0.141609_wp, 0.42158_wp, -0.145864_wp, -0.044507_wp, -0.019486_wp, &
                      -0.046969_wp, -0.604353_wp, 0.023802_wp, 0.038633_wp, -0.04461_wp, &
                      0.039344_wp, 0.025404_wp, 0.021645_wp, 0.016655_wp, 0.017719_wp, &
                      0.006981_wp, 0.016512_wp, 0.022736_wp, 0.126384_wp, 0.134998_wp, &
                      0.135006_wp]
         case ("mbis")
            charges = [ &
                      -0.229643_wp, 0.011973_wp, -0.229641_wp, -0.207833_wp, -0.204766_wp, &
                      -0.207834_wp, -0.468726_wp, 0.105845_wp, 0.105327_wp, 0.091165_wp, &
                      0.105323_wp, 0.105847_wp, 0.104979_wp, 0.101337_wp, 0.102173_wp, &
                      0.10388_wp, 0.10134_wp, 0.104977_wp, 0.133286_wp, 0.135496_wp, &
                      0.135497_wp]
         case ("chelpg")
            charges = [ &
                      -0.154981_wp, 0.439985_wp, -0.159338_wp, -0.0413_wp, -0.024217_wp, &
                      -0.043848_wp, -0.626062_wp, 0.026668_wp, 0.041829_wp, -0.046528_wp, &
                      0.042551_wp, 0.028307_wp, 0.022157_wp, 0.016628_wp, 0.01921_wp, &
                      0.007733_wp, 0.016501_wp, 0.023269_wp, 0.131512_wp, 0.139964_wp, &
                      0.139963_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.040183_wp, -0.005886_wp, -0.040182_wp, -0.038974_wp, -0.039202_wp, &
                      -0.038976_wp, -0.073257_wp, 0.020236_wp, 0.017483_wp, 0.017466_wp, &
                      0.017481_wp, 0.020237_wp, 0.021296_wp, 0.018096_wp, 0.018152_wp, &
                      0.021269_wp, 0.018097_wp, 0.021295_wp, 0.020816_wp, 0.022369_wp, &
                      0.022369_wp]
         case ("resp")
            charges = [ &
                      -0.149077_wp, 0.44212_wp, -0.153301_wp, -0.044928_wp, -0.020441_wp, &
                      -0.04741_wp, -0.62779_wp, 0.024598_wp, 0.04071_wp, -0.04643_wp, &
                      0.041415_wp, 0.026192_wp, 0.021643_wp, 0.017345_wp, 0.0187_wp, &
                      0.006649_wp, 0.017204_wp, 0.022749_wp, 0.130614_wp, 0.139717_wp, &
                      0.13972_wp]
         case ("mbis")
            charges = [ &
                      -0.232762_wp, 0.015801_wp, -0.23276_wp, -0.209906_wp, -0.206758_wp, &
                      -0.209908_wp, -0.477655_wp, 0.10689_wp, 0.106603_wp, 0.09209_wp, &
                      0.106599_wp, 0.106892_wp, 0.105747_wp, 0.102549_wp, 0.103355_wp, &
                      0.104641_wp, 0.102553_wp, 0.105745_wp, 0.135188_wp, 0.137548_wp, &
                      0.137549_wp]
         case ("chelpg")
            charges = [ &
                      -0.162706_wp, 0.4607_wp, -0.167023_wp, -0.041602_wp, -0.025325_wp, &
                      -0.044166_wp, -0.649626_wp, 0.027522_wp, 0.043973_wp, -0.048352_wp, &
                      0.044686_wp, 0.029151_wp, 0.022149_wp, 0.01731_wp, 0.020233_wp, &
                      0.007425_wp, 0.017185_wp, 0.023275_wp, 0.135776_wp, 0.144711_wp, &
                      0.144704_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.040699_wp, -0.005684_wp, -0.040698_wp, -0.039577_wp, -0.039814_wp, &
                      -0.039579_wp, -0.0747_wp, 0.020467_wp, 0.017869_wp, 0.018303_wp, &
                      0.017867_wp, 0.020468_wp, 0.021357_wp, 0.018565_wp, 0.018549_wp, &
                      0.021356_wp, 0.018566_wp, 0.021356_wp, 0.02091_wp, 0.022559_wp, &
                      0.022559_wp]
         case ("resp")
            charges = [ &
                      -0.157009_wp, 0.464775_wp, -0.161194_wp, -0.045079_wp, -0.021335_wp, &
                      -0.047582_wp, -0.652691_wp, 0.025286_wp, 0.042866_wp, -0.04834_wp, &
                      0.043564_wp, 0.026871_wp, 0.021383_wp, 0.018095_wp, 0.01971_wp, &
                      0.00608_wp, 0.017956_wp, 0.022504_wp, 0.134924_wp, 0.14461_wp, &
                      0.144607_wp]
         case ("mbis")
            charges = [ &
                      -0.235931_wp, 0.019815_wp, -0.235929_wp, -0.211983_wp, -0.208737_wp, &
                      -0.211986_wp, -0.486847_wp, 0.107891_wp, 0.107936_wp, 0.09329_wp, &
                      0.107932_wp, 0.107892_wp, 0.106364_wp, 0.103914_wp, 0.104599_wp, &
                      0.105273_wp, 0.103918_wp, 0.106362_wp, 0.137041_wp, 0.139594_wp, &
                      0.139594_wp]
         case ("chelpg")
            charges = [ &
                      -0.170891_wp, 0.483547_wp, -0.175162_wp, -0.041612_wp, -0.026377_wp, &
                      -0.044194_wp, -0.674627_wp, 0.028263_wp, 0.04619_wp, -0.050276_wp, &
                      0.046894_wp, 0.029881_wp, 0.021875_wp, 0.018044_wp, 0.021282_wp, &
                      0.006879_wp, 0.017921_wp, 0.023017_wp, 0.140111_wp, 0.149623_wp, &
                      0.14961_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (143)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.331538_wp, 0.156848_wp, 0.037841_wp, -0.105285_wp, -0.002638_wp, &
                      0.128773_wp, 0.033856_wp, 0.03582_wp, 0.046321_wp]
         case ("resp")
            charges = [ &
                      -0.48134_wp, 0.32963_wp, 0.042506_wp, -0.306697_wp, -0.088338_wp, &
                      0.255593_wp, 0.072121_wp, 0.063193_wp, 0.113334_wp]
         case ("mbis")
            charges = [ &
                      -0.572127_wp, 0.551064_wp, 0.024745_wp, -0.554225_wp, -0.15559_wp, &
                      0.348644_wp, 0.103445_wp, 0.105331_wp, 0.14871_wp]
         case ("chelpg")
            charges = [ &
                      -0.482801_wp, 0.334832_wp, 0.040851_wp, -0.309655_wp, -0.093577_wp, &
                      0.256899_wp, 0.073717_wp, 0.06483_wp, 0.114904_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.404289_wp, 0.157347_wp, 0.055208_wp, -0.088758_wp, -0.001224_wp, &
                      0.156759_wp, 0.040909_wp, 0.047549_wp, 0.036496_wp]
         case ("resp")
            charges = [ &
                      -0.600485_wp, 0.358038_wp, 0.069287_wp, -0.306417_wp, -0.088892_wp, &
                      0.303678_wp, 0.083546_wp, 0.084009_wp, 0.097236_wp]
         case ("mbis")
            charges = [ &
                      -0.674795_wp, 0.564411_wp, 0.056039_wp, -0.551649_wp, -0.16386_wp, &
                      0.394639_wp, 0.116033_wp, 0.126412_wp, 0.132768_wp]
         case ("chelpg")
            charges = [ &
                      -0.601903_wp, 0.363122_wp, 0.067664_wp, -0.309245_wp, -0.094348_wp, &
                      0.304956_wp, 0.085197_wp, 0.085701_wp, 0.098856_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.404802_wp, 0.157341_wp, 0.055326_wp, -0.088622_wp, -0.001209_wp, &
                      0.156959_wp, 0.040961_wp, 0.04764_wp, 0.036404_wp]
         case ("resp")
            charges = [ &
                      -0.601333_wp, 0.358238_wp, 0.069466_wp, -0.306378_wp, -0.088857_wp, &
                      0.304012_wp, 0.083622_wp, 0.084169_wp, 0.097062_wp]
         case ("mbis")
            charges = [ &
                      -0.675519_wp, 0.564485_wp, 0.056251_wp, -0.551573_wp, -0.163897_wp, &
                      0.394949_wp, 0.116121_wp, 0.126571_wp, 0.13261_wp]
         case ("chelpg")
            charges = [ &
                      -0.60276_wp, 0.363323_wp, 0.067846_wp, -0.309201_wp, -0.094317_wp, &
                      0.30529_wp, 0.085274_wp, 0.085863_wp, 0.098682_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (144)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.014043_wp, -0.023489_wp, -0.016725_wp, 0.037584_wp, -0.016732_wp, &
                      -0.023492_wp, 0.046277_wp, 0.046355_wp, 0.050659_wp, 0.234916_wp, &
                      0.050658_wp, 0.046353_wp, -0.209097_wp, -0.209225_wp]
         case ("resp")
            charges = [ &
                      -0.157697_wp, -0.06663_wp, -0.214139_wp, 0.145073_wp, -0.213797_wp, &
                      -0.067038_wp, 0.144704_wp, 0.131064_wp, 0.169392_wp, 0.626891_wp, &
                      0.169174_wp, 0.131414_wp, -0.398987_wp, -0.399424_wp]
         case ("mbis")
            charges = [ &
                      -0.102784_wp, -0.132166_wp, -0.14159_wp, 0.021789_wp, -0.141612_wp, &
                      -0.132165_wp, 0.143452_wp, 0.148067_wp, 0.174744_wp, 0.642674_wp, &
                      0.174728_wp, 0.148064_wp, -0.401523_wp, -0.401681_wp]
         case ("chelpg")
            charges = [ &
                      -0.17254_wp, -0.051242_wp, -0.232063_wp, 0.162818_wp, -0.231665_wp, &
                      -0.051687_wp, 0.147343_wp, 0.128519_wp, 0.173021_wp, 0.624676_wp, &
                      0.17278_wp, 0.128878_wp, -0.399198_wp, -0.399638_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.00837_wp, -0.021507_wp, -0.017418_wp, 0.032827_wp, -0.017417_wp, &
                      -0.021508_wp, 0.058459_wp, 0.057328_wp, 0.054381_wp, 0.239756_wp, &
                      0.05439_wp, 0.057326_wp, -0.234043_wp, -0.234206_wp]
         case ("resp")
            charges = [ &
                      -0.159032_wp, -0.070239_wp, -0.245808_wp, 0.173301_wp, -0.245334_wp, &
                      -0.070654_wp, 0.165838_wp, 0.152025_wp, 0.190747_wp, 0.664494_wp, &
                      0.190472_wp, 0.152363_wp, -0.448844_wp, -0.449329_wp]
         case ("mbis")
            charges = [ &
                      -0.100328_wp, -0.142213_wp, -0.142724_wp, 0.004457_wp, -0.142734_wp, &
                      -0.142212_wp, 0.163728_wp, 0.168496_wp, 0.182975_wp, 0.682438_wp, &
                      0.182969_wp, 0.168494_wp, -0.441573_wp, -0.441776_wp]
         case ("chelpg")
            charges = [ &
                      -0.174068_wp, -0.054633_wp, -0.264124_wp, 0.191558_wp, -0.263593_wp, &
                      -0.055085_wp, 0.168505_wp, 0.149457_wp, 0.194455_wp, 0.662099_wp, &
                      0.194159_wp, 0.149805_wp, -0.449025_wp, -0.449512_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.008031_wp, -0.021362_wp, -0.017466_wp, 0.032586_wp, -0.017465_wp, &
                      -0.021362_wp, 0.058944_wp, 0.057738_wp, 0.054471_wp, 0.239851_wp, &
                      0.054479_wp, 0.057737_wp, -0.23498_wp, -0.235144_wp]
         case ("resp")
            charges = [ &
                      -0.158767_wp, -0.070232_wp, -0.247039_wp, 0.174276_wp, -0.246561_wp, &
                      -0.070646_wp, 0.166581_wp, 0.152755_wp, 0.191465_wp, 0.665637_wp, &
                      0.191188_wp, 0.153093_wp, -0.450632_wp, -0.451118_wp]
         case ("mbis")
            charges = [ &
                      -0.099987_wp, -0.142456_wp, -0.142779_wp, 0.003698_wp, -0.142788_wp, &
                      -0.142455_wp, 0.16446_wp, 0.169215_wp, 0.183197_wp, 0.683752_wp, &
                      0.183192_wp, 0.169213_wp, -0.44303_wp, -0.443236_wp]
         case ("chelpg")
            charges = [ &
                      -0.173821_wp, -0.054609_wp, -0.265378_wp, 0.192563_wp, -0.264843_wp, &
                      -0.055061_wp, 0.169252_wp, 0.150187_wp, 0.195179_wp, 0.663233_wp, &
                      0.19488_wp, 0.150534_wp, -0.450815_wp, -0.451303_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (145)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.033692_wp, 0.055121_wp, -0.06412_wp, 0.054628_wp, 0.236057_wp, &
                      -0.224367_wp, -0.206338_wp, 0.039258_wp, 0.038059_wp, 0.038008_wp]
         case ("resp")
            charges = [ &
                      -0.008501_wp, 0.071142_wp, -0.247415_wp, 0.06961_wp, 0.71391_wp, &
                      -0.455314_wp, -0.411015_wp, 0.092553_wp, 0.088096_wp, 0.086935_wp]
         case ("mbis")
            charges = [ &
                      -0.195278_wp, 0.158458_wp, -0.444743_wp, 0.158616_wp, 0.649856_wp, &
                      -0.418433_wp, -0.388011_wp, 0.154829_wp, 0.162054_wp, 0.162649_wp]
         case ("chelpg")
            charges = [ &
                      -0.008729_wp, 0.071727_wp, -0.258255_wp, 0.070229_wp, 0.716668_wp, &
                      -0.456477_wp, -0.411623_wp, 0.095624_wp, 0.090984_wp, 0.089852_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.043193_wp, 0.068194_wp, -0.062414_wp, 0.067142_wp, 0.245456_wp, &
                      -0.254248_wp, -0.231454_wp, 0.049452_wp, 0.037563_wp, 0.037114_wp]
         case ("resp")
            charges = [ &
                      0.041347_wp, 0.079935_wp, -0.263806_wp, 0.077796_wp, 0.759038_wp, &
                      -0.513683_wp, -0.464074_wp, 0.111279_wp, 0.086791_wp, 0.085376_wp]
         case ("mbis")
            charges = [ &
                      -0.202743_wp, 0.179721_wp, -0.447913_wp, 0.178953_wp, 0.697324_wp, &
                      -0.467421_wp, -0.430725_wp, 0.172186_wp, 0.160362_wp, 0.160254_wp]
         case ("chelpg")
            charges = [ &
                      0.04398_wp, 0.079823_wp, -0.276557_wp, 0.077711_wp, 0.760858_wp, &
                      -0.514659_wp, -0.464528_wp, 0.114683_wp, 0.090039_wp, 0.088651_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.043625_wp, 0.068794_wp, -0.062322_wp, 0.067716_wp, 0.245818_wp, &
                      -0.255525_wp, -0.232593_wp, 0.04996_wp, 0.037495_wp, 0.037029_wp]
         case ("resp")
            charges = [ &
                      0.043782_wp, 0.080338_wp, -0.264585_wp, 0.078169_wp, 0.76069_wp, &
                      -0.516126_wp, -0.466402_wp, 0.112236_wp, 0.086661_wp, 0.085237_wp]
         case ("mbis")
            charges = [ &
                      -0.203058_wp, 0.180689_wp, -0.447931_wp, 0.17988_wp, 0.699277_wp, &
                      -0.469499_wp, -0.432625_wp, 0.173036_wp, 0.160184_wp, 0.160045_wp]
         case ("chelpg")
            charges = [ &
                      0.046533_wp, 0.080198_wp, -0.277416_wp, 0.078056_wp, 0.762471_wp, &
                      -0.517093_wp, -0.466853_wp, 0.115655_wp, 0.089923_wp, 0.088526_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (146)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.009301_wp, 0.059973_wp, 0.065283_wp, 0.059977_wp, 0.239757_wp, &
                      -0.217087_wp, -0.217204_wp]
         case ("resp")
            charges = [ &
                      -0.462814_wp, 0.1778_wp, 0.184304_wp, 0.178086_wp, 0.778043_wp, &
                      -0.427937_wp, -0.427482_wp]
         case ("mbis")
            charges = [ &
                      -0.440739_wp, 0.189676_wp, 0.189371_wp, 0.189671_wp, 0.691213_wp, &
                      -0.409525_wp, -0.409668_wp]
         case ("chelpg")
            charges = [ &
                      -0.477849_wp, 0.181727_wp, 0.188011_wp, 0.182028_wp, 0.783982_wp, &
                      -0.429186_wp, -0.428713_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.019982_wp, 0.070639_wp, 0.080127_wp, 0.070654_wp, 0.249661_wp, &
                      -0.245455_wp, -0.245609_wp]
         case ("resp")
            charges = [ &
                      -0.427332_wp, 0.186779_wp, 0.197509_wp, 0.187093_wp, 0.819107_wp, &
                      -0.481839_wp, -0.481317_wp]
         case ("mbis")
            charges = [ &
                      -0.453643_wp, 0.206774_wp, 0.213213_wp, 0.206787_wp, 0.739909_wp, &
                      -0.456426_wp, -0.456615_wp]
         case ("chelpg")
            charges = [ &
                      -0.442407_wp, 0.190716_wp, 0.201227_wp, 0.191044_wp, 0.825062_wp, &
                      -0.48309_wp, -0.482551_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.020354_wp, 0.071_wp, 0.080635_wp, 0.071016_wp, 0.249959_wp, &
                      -0.246405_wp, -0.24656_wp]
         case ("resp")
            charges = [ &
                      -0.425837_wp, 0.187025_wp, 0.19792_wp, 0.18734_wp, 0.82023_wp, &
                      -0.483601_wp, -0.483077_wp]
         case ("mbis")
            charges = [ &
                      -0.454018_wp, 0.207343_wp, 0.214019_wp, 0.207357_wp, 0.74144_wp, &
                      -0.457975_wp, -0.458167_wp]
         case ("chelpg")
            charges = [ &
                      -0.44091_wp, 0.190962_wp, 0.201638_wp, 0.191292_wp, 0.826183_wp, &
                      -0.484853_wp, -0.484311_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (147)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.06815_wp, 0.026285_wp, -0.032091_wp, -0.017479_wp, -0.029661_wp, &
                      -0.020488_wp, 0.032766_wp, 0.23639_wp, -0.199249_wp, -0.211985_wp, &
                      0.035056_wp, 0.038933_wp, 0.031954_wp, 0.040115_wp, 0.044694_wp, &
                      0.044448_wp, 0.048457_wp]
         case ("resp")
            charges = [ &
                      -0.565724_wp, 0.262362_wp, -0.21092_wp, -0.127681_wp, -0.110917_wp, &
                      -0.203257_wp, -0.054893_wp, 0.729239_wp, -0.419894_wp, -0.429712_wp, &
                      0.184701_wp, 0.180926_wp, 0.147303_wp, 0.150118_wp, 0.144573_wp, &
                      0.138516_wp, 0.185259_wp]
         case ("mbis")
            charges = [ &
                      -0.479303_wp, 0.147223_wp, -0.194853_wp, -0.093142_wp, -0.147773_wp, &
                      -0.132119_wp, -0.049076_wp, 0.668523_wp, -0.397366_wp, -0.411473_wp, &
                      0.16697_wp, 0.1717_wp, 0.147282_wp, 0.146125_wp, 0.140349_wp, &
                      0.147342_wp, 0.169588_wp]
         case ("chelpg")
            charges = [ &
                      -0.591949_wp, 0.285049_wp, -0.223435_wp, -0.118974_wp, -0.119297_wp, &
                      -0.193704_wp, -0.071399_wp, 0.735884_wp, -0.421571_wp, -0.431005_wp, &
                      0.190605_wp, 0.186904_wp, 0.153264_wp, 0.152231_wp, 0.143403_wp, &
                      0.139908_wp, 0.184085_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068905_wp, 0.026555_wp, -0.030947_wp, -0.01394_wp, -0.029628_wp, &
                      -0.023128_wp, 0.028132_wp, 0.241829_wp, -0.219152_wp, -0.237918_wp, &
                      0.035853_wp, 0.037359_wp, 0.040637_wp, 0.05113_wp, 0.055755_wp, &
                      0.054114_wp, 0.052248_wp]
         case ("resp")
            charges = [ &
                      -0.592401_wp, 0.27226_wp, -0.230169_wp, -0.129715_wp, -0.119775_wp, &
                      -0.233106_wp, -0.04617_wp, 0.782001_wp, -0.468809_wp, -0.484502_wp, &
                      0.192232_wp, 0.185295_wp, 0.167492_wp, 0.175156_wp, 0.165117_wp, &
                      0.158962_wp, 0.206132_wp]
         case ("mbis")
            charges = [ &
                      -0.492488_wp, 0.15599_wp, -0.208586_wp, -0.093498_wp, -0.1591_wp, &
                      -0.138772_wp, -0.065876_wp, 0.709958_wp, -0.43209_wp, -0.453724_wp, &
                      0.170103_wp, 0.170463_wp, 0.164618_wp, 0.167614_wp, 0.15972_wp, &
                      0.1664_wp, 0.179266_wp]
         case ("chelpg")
            charges = [ &
                      -0.618674_wp, 0.294651_wp, -0.242501_wp, -0.121189_wp, -0.127938_wp, &
                      -0.22398_wp, -0.062151_wp, 0.788448_wp, -0.470448_wp, -0.485762_wp, &
                      0.198156_wp, 0.191291_wp, 0.173485_wp, 0.177251_wp, 0.163974_wp, &
                      0.160335_wp, 0.205052_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068932_wp, 0.026557_wp, -0.030775_wp, -0.01361_wp, -0.029547_wp, &
                      -0.023305_wp, 0.027827_wp, 0.241996_wp, -0.220176_wp, -0.239237_wp, &
                      0.035819_wp, 0.03723_wp, 0.041117_wp, 0.051723_wp, 0.056354_wp, &
                      0.0546_wp, 0.052357_wp]
         case ("resp")
            charges = [ &
                      -0.593415_wp, 0.272621_wp, -0.230919_wp, -0.129468_wp, -0.120057_wp, &
                      -0.234674_wp, -0.045809_wp, 0.784305_wp, -0.471208_wp, -0.487165_wp, &
                      0.192406_wp, 0.185359_wp, 0.168534_wp, 0.176402_wp, 0.166109_wp, &
                      0.159932_wp, 0.207048_wp]
         case ("mbis")
            charges = [ &
                      -0.492963_wp, 0.15633_wp, -0.20905_wp, -0.093247_wp, -0.159504_wp, &
                      -0.139158_wp, -0.066832_wp, 0.711839_wp, -0.433819_wp, -0.455815_wp, &
                      0.170102_wp, 0.170284_wp, 0.16554_wp, 0.168687_wp, 0.160677_wp, &
                      0.167305_wp, 0.179621_wp]
         case ("chelpg")
            charges = [ &
                      -0.619689_wp, 0.295_wp, -0.243247_wp, -0.120946_wp, -0.128212_wp, &
                      -0.225566_wp, -0.061767_wp, 0.790743_wp, -0.472845_wp, -0.488424_wp, &
                      0.19833_wp, 0.191355_wp, 0.174528_wp, 0.178497_wp, 0.164966_wp, &
                      0.161304_wp, 0.205972_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (148)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070368_wp, -0.034194_wp, -0.036417_wp, -0.036256_wp, -0.036168_wp, &
                      -0.036256_wp, -0.036417_wp, -0.034194_wp, -0.070368_wp, 0.021429_wp, &
                      0.02143_wp, 0.023255_wp, 0.019084_wp, 0.019094_wp, 0.018209_wp, &
                      0.018201_wp, 0.0182_wp, 0.01822_wp, 0.018219_wp, 0.018188_wp, &
                      0.018197_wp, 0.018224_wp, 0.018218_wp, 0.018192_wp, 0.019082_wp, &
                      0.019095_wp, 0.021435_wp, 0.021424_wp, 0.023255_wp]
         case ("resp")
            charges = [ &
                      -0.350725_wp, 0.223889_wp, -0.011659_wp, -0.07202_wp, 0.14959_wp, &
                      -0.071889_wp, -0.014297_wp, 0.228169_wp, -0.344526_wp, 0.075893_wp, &
                      0.07579_wp, 0.079588_wp, -0.036538_wp, -0.036443_wp, 0.000494_wp, &
                      1e-05_wp, 0.004663_wp, 0.004917_wp, -0.032387_wp, -0.032843_wp, &
                      0.005129_wp, 0.005021_wp, 0.000555_wp, 0.000519_wp, -0.038209_wp, &
                      -0.038438_wp, 0.07419_wp, 0.073651_wp, 0.077908_wp]
         case ("mbis")
            charges = [ &
                      -0.421681_wp, -0.160074_wp, -0.211327_wp, -0.213999_wp, -0.207118_wp, &
                      -0.213997_wp, -0.211327_wp, -0.160075_wp, -0.421681_wp, 0.127481_wp, &
                      0.127473_wp, 0.129659_wp, 0.099798_wp, 0.099857_wp, 0.105605_wp, &
                      0.10555_wp, 0.105279_wp, 0.105306_wp, 0.104673_wp, 0.104609_wp, &
                      0.105277_wp, 0.105308_wp, 0.105599_wp, 0.105556_wp, 0.099812_wp, &
                      0.099841_wp, 0.127493_wp, 0.12746_wp, 0.129658_wp]
         case ("chelpg")
            charges = [ &
                      -0.37096_wp, 0.240402_wp, -0.013779_wp, -0.087292_wp, 0.173801_wp, &
                      -0.086909_wp, -0.016888_wp, 0.244931_wp, -0.36476_wp, 0.080386_wp, &
                      0.080278_wp, 0.084056_wp, -0.039422_wp, -0.0393_wp, 0.000633_wp, &
                      0.000168_wp, 0.007319_wp, 0.007559_wp, -0.036695_wp, -0.037185_wp, &
                      0.007769_wp, 0.007636_wp, 0.000799_wp, 0.000772_wp, -0.041145_wp, &
                      -0.041348_wp, 0.078666_wp, 0.078147_wp, 0.082361_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071695_wp, -0.034661_wp, -0.03688_wp, -0.036749_wp, -0.036681_wp, &
                      -0.036749_wp, -0.03688_wp, -0.034661_wp, -0.071695_wp, 0.021691_wp, &
                      0.021694_wp, 0.02349_wp, 0.019435_wp, 0.019441_wp, 0.018542_wp, &
                      0.018537_wp, 0.018506_wp, 0.018521_wp, 0.018488_wp, 0.018464_wp, &
                      0.018502_wp, 0.018525_wp, 0.018549_wp, 0.01853_wp, 0.019433_wp, &
                      0.019443_wp, 0.021696_wp, 0.021688_wp, 0.02349_wp]
         case ("resp")
            charges = [ &
                      -0.36151_wp, 0.234106_wp, -0.009687_wp, -0.07144_wp, 0.157452_wp, &
                      -0.071629_wp, -0.01227_wp, 0.238254_wp, -0.355172_wp, 0.077794_wp, &
                      0.077699_wp, 0.081572_wp, -0.038762_wp, -0.038669_wp, -0.000625_wp, &
                      -0.001112_wp, 0.003449_wp, 0.003691_wp, -0.034924_wp, -0.035394_wp, &
                      0.004003_wp, 0.003891_wp, -0.000554_wp, -0.000578_wp, -0.040408_wp, &
                      -0.040634_wp, 0.076064_wp, 0.075527_wp, 0.079866_wp]
         case ("mbis")
            charges = [ &
                      -0.428831_wp, -0.160364_wp, -0.213284_wp, -0.215972_wp, -0.20899_wp, &
                      -0.21597_wp, -0.213285_wp, -0.160366_wp, -0.428831_wp, 0.129308_wp, &
                      0.129303_wp, 0.131447_wp, 0.100697_wp, 0.100747_wp, 0.10666_wp, &
                      0.106612_wp, 0.106287_wp, 0.106306_wp, 0.105614_wp, 0.105564_wp, &
                      0.106285_wp, 0.106308_wp, 0.106651_wp, 0.106622_wp, 0.10071_wp, &
                      0.100733_wp, 0.129318_wp, 0.129291_wp, 0.131447_wp]
         case ("chelpg")
            charges = [ &
                      -0.381878_wp, 0.250713_wp, -0.011669_wp, -0.086839_wp, 0.181943_wp, &
                      -0.086807_wp, -0.014703_wp, 0.255099_wp, -0.37554_wp, 0.082314_wp, &
                      0.082214_wp, 0.086074_wp, -0.041672_wp, -0.041551_wp, -0.00053_wp, &
                      -0.000998_wp, 0.00611_wp, 0.00634_wp, -0.039298_wp, -0.039803_wp, &
                      0.006656_wp, 0.006518_wp, -0.000357_wp, -0.000372_wp, -0.043368_wp, &
                      -0.043567_wp, 0.080567_wp, 0.080051_wp, 0.084353_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073156_wp, -0.035136_wp, -0.037364_wp, -0.037266_wp, -0.037223_wp, &
                      -0.037266_wp, -0.037364_wp, -0.035137_wp, -0.073155_wp, 0.021946_wp, &
                      0.021951_wp, 0.023687_wp, 0.019835_wp, 0.019834_wp, 0.018906_wp, &
                      0.018907_wp, 0.018851_wp, 0.018859_wp, 0.01877_wp, 0.018758_wp, &
                      0.018848_wp, 0.018862_wp, 0.018911_wp, 0.018903_wp, 0.019832_wp, &
                      0.019836_wp, 0.021949_wp, 0.021947_wp, 0.023687_wp]
         case ("resp")
            charges = [ &
                      -0.373376_wp, 0.24618_wp, -0.007017_wp, -0.070253_wp, 0.166805_wp, &
                      -0.070819_wp, -0.009534_wp, 0.250175_wp, -0.366879_wp, 0.079775_wp, &
                      0.079689_wp, 0.083612_wp, -0.041455_wp, -0.041366_wp, -0.002115_wp, &
                      -0.002601_wp, 0.001851_wp, 0.002074_wp, -0.038046_wp, -0.038525_wp, &
                      0.002508_wp, 0.002385_wp, -0.002034_wp, -0.002038_wp, -0.043071_wp, &
                      -0.043298_wp, 0.078013_wp, 0.077482_wp, 0.081876_wp]
         case ("mbis")
            charges = [ &
                      -0.436609_wp, -0.160586_wp, -0.215362_wp, -0.218075_wp, -0.210976_wp, &
                      -0.218073_wp, -0.215362_wp, -0.160587_wp, -0.436608_wp, 0.131254_wp, &
                      0.131252_wp, 0.133298_wp, 0.101684_wp, 0.101721_wp, 0.107794_wp, &
                      0.107757_wp, 0.107388_wp, 0.107392_wp, 0.106602_wp, 0.106574_wp, &
                      0.107387_wp, 0.107394_wp, 0.107779_wp, 0.107772_wp, 0.101695_wp, &
                      0.101708_wp, 0.131261_wp, 0.131243_wp, 0.133298_wp]
         case ("chelpg")
            charges = [ &
                      -0.393878_wp, 0.262868_wp, -0.008803_wp, -0.085753_wp, 0.191572_wp, &
                      -0.086133_wp, -0.011747_wp, 0.267087_wp, -0.387382_wp, 0.084322_wp, &
                      0.084232_wp, 0.088147_wp, -0.044392_wp, -0.044275_wp, -0.002084_wp, &
                      -0.002551_wp, 0.004507_wp, 0.004717_wp, -0.042491_wp, -0.043006_wp, &
                      0.005164_wp, 0.005015_wp, -0.001904_wp, -0.001898_wp, -0.046055_wp, &
                      -0.046255_wp, 0.082544_wp, 0.082033_wp, 0.086399_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (149)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036221_wp, 0.018229_wp, 0.018237_wp, -0.036421_wp, 0.018202_wp, &
                      0.0182_wp, -0.034191_wp, 0.019093_wp, 0.019102_wp, -0.070364_wp, &
                      0.021435_wp, 0.023266_wp, 0.021438_wp, -0.036221_wp, 0.01823_wp, &
                      0.018237_wp, -0.036421_wp, 0.018201_wp, 0.018201_wp, -0.034191_wp, &
                      0.019094_wp, 0.019101_wp, -0.070364_wp, 0.021435_wp, 0.021437_wp, &
                      0.023266_wp]
         case ("resp")
            charges = [ &
                      0.006711_wp, -0.011475_wp, -0.011614_wp, 0.000559_wp, -0.000246_wp, &
                      -0.001195_wp, 0.200392_wp, -0.034359_wp, -0.033675_wp, -0.336273_wp, &
                      0.072131_wp, 0.079494_wp, 0.071426_wp, -0.001486_wp, -0.01055_wp, &
                      -0.009303_wp, 0.001487_wp, 0.000309_wp, 0.000636_wp, 0.188329_wp, &
                      -0.029862_wp, -0.030943_wp, -0.330721_wp, 0.069963_wp, 0.071643_wp, &
                      0.078622_wp]
         case ("mbis")
            charges = [ &
                      -0.212651_wp, 0.105133_wp, 0.105154_wp, -0.211285_wp, 0.105581_wp, &
                      0.105586_wp, -0.160212_wp, 0.099831_wp, 0.099881_wp, -0.4217_wp, &
                      0.127503_wp, 0.12968_wp, 0.127503_wp, -0.212652_wp, 0.105129_wp, &
                      0.105158_wp, -0.211283_wp, 0.105584_wp, 0.105583_wp, -0.160212_wp, &
                      0.099835_wp, 0.099877_wp, -0.421706_wp, 0.127506_wp, 0.127503_wp, &
                      0.129682_wp]
         case ("chelpg")
            charges = [ &
                      0.010558_wp, -0.01231_wp, -0.01244_wp, -0.007021_wp, 0.000863_wp, &
                      -8.4e-05_wp, 0.214349_wp, -0.037038_wp, -0.036345_wp, -0.351839_wp, &
                      0.075677_wp, 0.08292_wp, 0.074958_wp, 0.000362_wp, -0.010996_wp, &
                      -0.009746_wp, -0.005323_wp, 0.00131_wp, 0.001641_wp, 0.202247_wp, &
                      -0.032508_wp, -0.033572_wp, -0.346614_wp, 0.073576_wp, 0.07523_wp, &
                      0.082145_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036702_wp, 0.018531_wp, 0.018535_wp, -0.036875_wp, 0.018537_wp, &
                      0.018534_wp, -0.034652_wp, 0.019439_wp, 0.019444_wp, -0.071675_wp, &
                      0.021695_wp, 0.023496_wp, 0.021698_wp, -0.036702_wp, 0.018531_wp, &
                      0.018535_wp, -0.036875_wp, 0.018536_wp, 0.018535_wp, -0.034652_wp, &
                      0.01944_wp, 0.019443_wp, -0.071675_wp, 0.021695_wp, 0.021698_wp, &
                      0.023496_wp]
         case ("resp")
            charges = [ &
                      0.009606_wp, -0.013074_wp, -0.013217_wp, 0.002775_wp, -0.001283_wp, &
                      -0.002272_wp, 0.20949_wp, -0.036376_wp, -0.035696_wp, -0.346506_wp, &
                      0.073912_wp, 0.081447_wp, 0.073201_wp, 0.000869_wp, -0.012071_wp, &
                      -0.010802_wp, 0.003852_wp, -0.000787_wp, -0.000472_wp, 0.197558_wp, &
                      -0.03193_wp, -0.033029_wp, -0.340768_wp, 0.071679_wp, 0.073377_wp, &
                      0.080518_wp]
         case ("mbis")
            charges = [ &
                      -0.214591_wp, 0.106124_wp, 0.106139_wp, -0.213216_wp, 0.106635_wp, &
                      0.106636_wp, -0.160493_wp, 0.100717_wp, 0.100758_wp, -0.428776_wp, &
                      0.129312_wp, 0.131445_wp, 0.129315_wp, -0.214591_wp, 0.10612_wp, &
                      0.106143_wp, -0.213215_wp, 0.106638_wp, 0.106634_wp, -0.160493_wp, &
                      0.10072_wp, 0.100754_wp, -0.428782_wp, 0.129315_wp, 0.129314_wp, &
                      0.131448_wp]
         case ("chelpg")
            charges = [ &
                      0.013605_wp, -0.013969_wp, -0.014103_wp, -0.004741_wp, -0.00021_wp, &
                      -0.001195_wp, 0.223514_wp, -0.039076_wp, -0.038388_wp, -0.362174_wp, &
                      0.077481_wp, 0.084899_wp, 0.076755_wp, 0.002761_wp, -0.012558_wp, &
                      -0.011285_wp, -0.002851_wp, 0.000171_wp, 0.000491_wp, 0.211557_wp, &
                      -0.0346_wp, -0.035681_wp, -0.356787_wp, 0.075319_wp, 0.076991_wp, &
                      0.084073_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037218_wp, 0.018876_wp, 0.018876_wp, -0.037356_wp, 0.018918_wp, &
                      0.018912_wp, -0.03513_wp, 0.019839_wp, 0.019838_wp, -0.073152_wp, &
                      0.021954_wp, 0.02369_wp, 0.021959_wp, -0.037218_wp, 0.018877_wp, &
                      0.018876_wp, -0.037356_wp, 0.018916_wp, 0.018914_wp, -0.035131_wp, &
                      0.019839_wp, 0.019837_wp, -0.073152_wp, 0.021954_wp, 0.021958_wp, &
                      0.02369_wp]
         case ("resp")
            charges = [ &
                      0.013438_wp, -0.015133_wp, -0.015284_wp, 0.005759_wp, -0.002693_wp, &
                      -0.003728_wp, 0.220536_wp, -0.038892_wp, -0.038219_wp, -0.357993_wp, &
                      0.075801_wp, 0.083492_wp, 0.075081_wp, 0.00405_wp, -0.014036_wp, &
                      -0.012743_wp, 0.007006_wp, -0.002264_wp, -0.001962_wp, 0.208761_wp, &
                      -0.034506_wp, -0.035628_wp, -0.352035_wp, 0.07349_wp, 0.07521_wp, &
                      0.082494_wp]
         case ("mbis")
            charges = [ &
                      -0.216707_wp, 0.107223_wp, 0.107231_wp, -0.215308_wp, 0.1078_wp, &
                      0.107797_wp, -0.160705_wp, 0.101705_wp, 0.101734_wp, -0.43664_wp, &
                      0.131282_wp, 0.133308_wp, 0.131287_wp, -0.216708_wp, 0.107221_wp, &
                      0.107234_wp, -0.215307_wp, 0.107802_wp, 0.107795_wp, -0.160705_wp, &
                      0.101707_wp, 0.101731_wp, -0.436646_wp, 0.131284_wp, 0.131287_wp, &
                      0.13331_wp]
         case ("chelpg")
            charges = [ &
                      0.01763_wp, -0.016109_wp, -0.01625_wp, -0.001638_wp, -0.001672_wp, &
                      -0.002706_wp, 0.234615_wp, -0.041614_wp, -0.040932_wp, -0.373766_wp, &
                      0.079393_wp, 0.086971_wp, 0.078659_wp, 0.006014_wp, -0.01458_wp, &
                      -0.013283_wp, 0.000468_wp, -0.001368_wp, -0.001061_wp, 0.22283_wp, &
                      -0.037201_wp, -0.038305_wp, -0.368187_wp, 0.077158_wp, 0.078852_wp, &
                      0.086083_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (150)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070879_wp, -0.035445_wp, -0.037811_wp, -0.036485_wp, -0.03845_wp, &
                      -0.037186_wp, -0.036287_wp, -0.038273_wp, -0.036405_wp, -0.03662_wp, &
                      -0.03749_wp, -0.036611_wp, -0.036907_wp, -0.034674_wp, -0.071047_wp, &
                      0.021483_wp, 0.022872_wp, 0.021569_wp, 0.019354_wp, 0.019388_wp, &
                      0.018181_wp, 0.019651_wp, 0.018772_wp, 0.018749_wp, 0.019509_wp, &
                      0.018399_wp, 0.018452_wp, 0.019504_wp, 0.018853_wp, 0.018132_wp, &
                      0.019273_wp, 0.017728_wp, 0.019017_wp, 0.01829_wp, 0.017973_wp, &
                      0.019344_wp, 0.018305_wp, 0.018819_wp, 0.01867_wp, 0.01837_wp, &
                      0.0185_wp, 0.018497_wp, 0.019331_wp, 0.019359_wp, 0.021561_wp, &
                      0.023124_wp, 0.021559_wp]
         case ("resp")
            charges = [ &
                      -0.486379_wp, 0.285969_wp, 0.049027_wp, -0.285835_wp, 0.081907_wp, &
                      0.069241_wp, -0.334079_wp, 0.16453_wp, -0.168786_wp, -0.078871_wp, &
                      0.133881_wp, -0.151782_wp, 0.027742_wp, 0.242195_wp, -0.373869_wp, &
                      0.104405_wp, 0.10948_wp, 0.107455_wp, -0.035728_wp, -0.040071_wp, &
                      -0.021874_wp, 0.002419_wp, 0.069892_wp, 0.060997_wp, 0.008312_wp, &
                      -0.026398_wp, 0.018155_wp, 0.009122_wp, 0.073399_wp, 0.078729_wp, &
                      -0.010817_wp, 0.012115_wp, 0.049873_wp, 0.045168_wp, 0.018611_wp, &
                      0.029659_wp, -0.01013_wp, -0.022841_wp, 0.022678_wp, 0.02519_wp, &
                      -0.008405_wp, -0.004665_wp, -0.041346_wp, -0.041161_wp, 0.079017_wp, &
                      0.08405_wp, 0.079819_wp]
         case ("mbis")
            charges = [ &
                      -0.427242_wp, -0.15874_wp, -0.209156_wp, -0.225164_wp, -0.199711_wp, &
                      -0.208319_wp, -0.222224_wp, -0.193465_wp, -0.215857_wp, -0.212264_wp, &
                      -0.199799_wp, -0.219404_wp, -0.210405_wp, -0.159966_wp, -0.421562_wp, &
                      0.128175_wp, 0.129953_wp, 0.128644_wp, 0.099681_wp, 0.104203_wp, &
                      0.105535_wp, 0.105791_wp, 0.104804_wp, 0.10468_wp, 0.102947_wp, &
                      0.10504_wp, 0.107664_wp, 0.104938_wp, 0.105403_wp, 0.104954_wp, &
                      0.108529_wp, 0.09963_wp, 0.104572_wp, 0.104596_wp, 0.103348_wp, &
                      0.105359_wp, 0.102743_wp, 0.105544_wp, 0.106652_wp, 0.105308_wp, &
                      0.105567_wp, 0.105275_wp, 0.099699_wp, 0.099855_wp, 0.127517_wp, &
                      0.129221_wp, 0.127467_wp]
         case ("chelpg")
            charges = [ &
                      -0.504853_wp, 0.29799_wp, 0.056743_wp, -0.309078_wp, 0.090254_wp, &
                      0.08285_wp, -0.361916_wp, 0.183925_wp, -0.181547_wp, -0.090435_wp, &
                      0.154207_wp, -0.172559_wp, 0.033968_wp, 0.25588_wp, -0.394543_wp, &
                      0.108166_wp, 0.114108_wp, 0.111525_wp, -0.038104_wp, -0.042893_wp, &
                      -0.023987_wp, 0.001524_wp, 0.075257_wp, 0.066401_wp, 0.007108_wp, &
                      -0.027649_wp, 0.015929_wp, 0.007241_wp, 0.07932_wp, 0.084256_wp, &
                      -0.014876_wp, 0.009605_wp, 0.053898_wp, 0.048502_wp, 0.021364_wp, &
                      0.033256_wp, -0.013279_wp, -0.026045_wp, 0.026404_wp, 0.02879_wp, &
                      -0.009848_wp, -0.006159_wp, -0.043861_wp, -0.043621_wp, 0.083521_wp, &
                      0.088988_wp, 0.084276_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072297_wp, -0.035668_wp, -0.038239_wp, -0.037152_wp, -0.038815_wp, &
                      -0.037552_wp, -0.036978_wp, -0.03855_wp, -0.037083_wp, -0.037232_wp, &
                      -0.037903_wp, -0.037149_wp, -0.03735_wp, -0.035132_wp, -0.072419_wp, &
                      0.021656_wp, 0.023013_wp, 0.021799_wp, 0.019723_wp, 0.02012_wp, &
                      0.018455_wp, 0.019769_wp, 0.018862_wp, 0.018862_wp, 0.019637_wp, &
                      0.018883_wp, 0.019218_wp, 0.019604_wp, 0.018893_wp, 0.018199_wp, &
                      0.019827_wp, 0.018323_wp, 0.019101_wp, 0.018424_wp, 0.018211_wp, &
                      0.019429_wp, 0.018658_wp, 0.019318_wp, 0.018898_wp, 0.01864_wp, &
                      0.018851_wp, 0.018867_wp, 0.019678_wp, 0.01969_wp, 0.021819_wp, &
                      0.023288_wp, 0.021823_wp]
         case ("resp")
            charges = [ &
                      -0.503683_wp, 0.300013_wp, 0.052251_wp, -0.298533_wp, 0.08553_wp, &
                      0.074598_wp, -0.350152_wp, 0.174572_wp, -0.176783_wp, -0.080417_wp, &
                      0.141846_wp, -0.156617_wp, 0.032207_wp, 0.253528_wp, -0.385629_wp, &
                      0.107481_wp, 0.112642_wp, 0.110921_wp, -0.037919_wp, -0.041827_wp, &
                      -0.023401_wp, 0.001742_wp, 0.072559_wp, 0.063419_wp, 0.008284_wp, &
                      -0.027258_wp, 0.019208_wp, 0.008635_wp, 0.076102_wp, 0.081884_wp, &
                      -0.011943_wp, 0.012838_wp, 0.051573_wp, 0.046752_wp, 0.0185_wp, &
                      0.029815_wp, -0.011293_wp, -0.024487_wp, 0.02253_wp, 0.0253_wp, &
                      -0.01005_wp, -0.006058_wp, -0.043845_wp, -0.043736_wp, 0.080957_wp, &
                      0.086076_wp, 0.081868_wp]
         case ("mbis")
            charges = [ &
                      -0.435055_wp, -0.158151_wp, -0.210669_wp, -0.228009_wp, -0.200811_wp, &
                      -0.209707_wp, -0.225257_wp, -0.194089_wp, -0.218541_wp, -0.214336_wp, &
                      -0.201284_wp, -0.221711_wp, -0.212289_wp, -0.160135_wp, -0.428958_wp, &
                      0.129972_wp, 0.131698_wp, 0.130541_wp, 0.100407_wp, 0.105618_wp, &
                      0.106373_wp, 0.106417_wp, 0.105624_wp, 0.105524_wp, 0.103496_wp, &
                      0.10618_wp, 0.109458_wp, 0.10544_wp, 0.106137_wp, 0.105751_wp, &
                      0.109753_wp, 0.100942_wp, 0.105358_wp, 0.105469_wp, 0.104312_wp, &
                      0.106047_wp, 0.103758_wp, 0.106832_wp, 0.107573_wp, 0.106323_wp, &
                      0.106642_wp, 0.106386_wp, 0.100583_wp, 0.100706_wp, 0.129398_wp, &
                      0.130941_wp, 0.12936_wp]
         case ("chelpg")
            charges = [ &
                      -0.522266_wp, 0.312092_wp, 0.060268_wp, -0.32214_wp, 0.094007_wp, &
                      0.088532_wp, -0.378411_wp, 0.194262_wp, -0.18979_wp, -0.092101_wp, &
                      0.162557_wp, -0.177892_wp, 0.038954_wp, 0.26721_wp, -0.406468_wp, &
                      0.111258_wp, 0.117303_wp, 0.115011_wp, -0.040322_wp, -0.044682_wp, &
                      -0.025595_wp, 0.000783_wp, 0.077992_wp, 0.068897_wp, 0.00705_wp, &
                      -0.028534_wp, 0.016912_wp, 0.006685_wp, 0.082107_wp, 0.087487_wp, &
                      -0.016069_wp, 0.010283_wp, 0.05567_wp, 0.050148_wp, 0.021287_wp, &
                      0.033454_wp, -0.014515_wp, -0.027759_wp, 0.026326_wp, 0.028972_wp, &
                      -0.011616_wp, -0.007677_wp, -0.046373_wp, -0.046208_wp, 0.085491_wp, &
                      0.091066_wp, 0.086355_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073753_wp, -0.035817_wp, -0.038666_wp, -0.037848_wp, -0.039152_wp, &
                      -0.037889_wp, -0.0377_wp, -0.03875_wp, -0.037803_wp, -0.037854_wp, &
                      -0.038278_wp, -0.037681_wp, -0.03776_wp, -0.035573_wp, -0.073827_wp, &
                      0.02171_wp, 0.023102_wp, 0.022021_wp, 0.02011_wp, 0.021103_wp, &
                      0.018768_wp, 0.019771_wp, 0.018832_wp, 0.01893_wp, 0.019678_wp, &
                      0.01955_wp, 0.020209_wp, 0.019579_wp, 0.018816_wp, 0.018127_wp, &
                      0.02051_wp, 0.019139_wp, 0.01906_wp, 0.018463_wp, 0.018478_wp, &
                      0.019388_wp, 0.019024_wp, 0.019957_wp, 0.019117_wp, 0.018862_wp, &
                      0.01922_wp, 0.019295_wp, 0.020035_wp, 0.020012_wp, 0.022065_wp, &
                      0.023372_wp, 0.022065_wp]
         case ("resp")
            charges = [ &
                      -0.521785_wp, 0.315527_wp, 0.055887_wp, -0.311933_wp, 0.089442_wp, &
                      0.080883_wp, -0.367509_wp, 0.185888_wp, -0.185458_wp, -0.081627_wp, &
                      0.150723_wp, -0.161567_wp, 0.037524_wp, 0.265894_wp, -0.39768_wp, &
                      0.110406_wp, 0.11582_wp, 0.114495_wp, -0.040477_wp, -0.043409_wp, &
                      -0.025086_wp, 0.000702_wp, 0.075124_wp, 0.065876_wp, 0.008015_wp, &
                      -0.027905_wp, 0.020537_wp, 0.007714_wp, 0.078776_wp, 0.084982_wp, &
                      -0.013135_wp, 0.013843_wp, 0.053143_wp, 0.048254_wp, 0.018259_wp, &
                      0.029621_wp, -0.012709_wp, -0.02615_wp, 0.022211_wp, 0.02517_wp, &
                      -0.011999_wp, -0.007669_wp, -0.046649_wp, -0.046657_wp, 0.082858_wp, &
                      0.087968_wp, 0.083862_wp]
         case ("mbis")
            charges = [ &
                      -0.44293_wp, -0.157519_wp, -0.212114_wp, -0.230853_wp, -0.20191_wp, &
                      -0.211082_wp, -0.228229_wp, -0.194725_wp, -0.22123_wp, -0.216339_wp, &
                      -0.202794_wp, -0.223978_wp, -0.214136_wp, -0.160214_wp, -0.436448_wp, &
                      0.131579_wp, 0.133377_wp, 0.132451_wp, 0.101125_wp, 0.107475_wp, &
                      0.10727_wp, 0.10682_wp, 0.106245_wp, 0.106304_wp, 0.10387_wp, &
                      0.107653_wp, 0.11167_wp, 0.105689_wp, 0.106666_wp, 0.106302_wp, &
                      0.111177_wp, 0.102642_wp, 0.105935_wp, 0.106192_wp, 0.105327_wp, &
                      0.106503_wp, 0.104782_wp, 0.108368_wp, 0.108466_wp, 0.107239_wp, &
                      0.107728_wp, 0.107586_wp, 0.101472_wp, 0.101525_wp, 0.13129_wp, &
                      0.132548_wp, 0.131245_wp]
         case ("chelpg")
            charges = [ &
                      -0.54047_wp, 0.327656_wp, 0.064231_wp, -0.335907_wp, 0.09805_wp, &
                      0.095157_wp, -0.396189_wp, 0.205868_wp, -0.198721_wp, -0.093398_wp, &
                      0.171816_wp, -0.183364_wp, 0.044857_wp, 0.279549_wp, -0.418675_wp, &
                      0.114196_wp, 0.120513_wp, 0.118602_wp, -0.042908_wp, -0.046298_wp, &
                      -0.027369_wp, -0.00033_wp, 0.080626_wp, 0.071428_wp, 0.006751_wp, &
                      -0.029206_wp, 0.018167_wp, 0.005689_wp, 0.084866_wp, 0.090659_wp, &
                      -0.017328_wp, 0.011243_wp, 0.057312_wp, 0.051714_wp, 0.021072_wp, &
                      0.033295_wp, -0.016005_wp, -0.029491_wp, 0.02608_wp, 0.028915_wp, &
                      -0.013704_wp, -0.009429_wp, -0.049187_wp, -0.049139_wp, 0.087421_wp, &
                      0.09301_wp, 0.088376_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (151)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.157213_wp, -0.03674_wp, -0.033785_wp, -0.03345_wp, -0.068326_wp, &
                      0.030913_wp, -0.25772_wp, 0.038377_wp, 0.038297_wp, 0.026553_wp, &
                      0.026553_wp, 0.019383_wp, 0.019417_wp, 0.024217_wp, 0.024923_wp, &
                      0.024181_wp]
         case ("resp")
            charges = [ &
                      0.502372_wp, -0.226751_wp, 0.111212_wp, 0.101443_wp, -0.359789_wp, &
                      -0.024448_wp, -0.48886_wp, 0.053806_wp, 0.053744_wp, 0.010644_wp, &
                      0.011495_wp, -0.003969_wp, -0.003049_wp, 0.083522_wp, 0.095377_wp, &
                      0.083251_wp]
         case ("mbis")
            charges = [ &
                      0.458498_wp, -0.387619_wp, -0.174776_wp, -0.169915_wp, -0.421033_wp, &
                      0.020769_wp, -0.461353_wp, 0.148381_wp, 0.148218_wp, 0.121514_wp, &
                      0.121045_wp, 0.10119_wp, 0.101141_wp, 0.131155_wp, 0.131744_wp, &
                      0.131047_wp]
         case ("chelpg")
            charges = [ &
                      0.512061_wp, -0.246651_wp, 0.121076_wp, 0.107494_wp, -0.376213_wp, &
                      -0.026352_wp, -0.491122_wp, 0.057873_wp, 0.057819_wp, 0.009097_wp, &
                      0.009948_wp, -0.004863_wp, -0.003949_wp, 0.087147_wp, 0.099756_wp, &
                      0.08688_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.168969_wp, -0.032961_wp, -0.035755_wp, -0.032126_wp, -0.071045_wp, &
                      0.045721_wp, -0.306245_wp, 0.049765_wp, 0.049644_wp, 0.022826_wp, &
                      0.023114_wp, 0.022696_wp, 0.022745_wp, 0.023685_wp, 0.02526_wp, &
                      0.023712_wp]
         case ("resp")
            charges = [ &
                      0.54931_wp, -0.226896_wp, 0.112552_wp, 0.137145_wp, -0.390777_wp, &
                      -0.006416_wp, -0.577158_wp, 0.068333_wp, 0.068254_wp, 0.001169_wp, &
                      0.002421_wp, -0.007186_wp, -0.006203_wp, 0.087618_wp, 0.100355_wp, &
                      0.087481_wp]
         case ("mbis")
            charges = [ &
                      0.498303_wp, -0.410636_wp, -0.171683_wp, -0.167097_wp, -0.435716_wp, &
                      0.04159_wp, -0.537834_wp, 0.169428_wp, 0.169189_wp, 0.115239_wp, &
                      0.115286_wp, 0.106287_wp, 0.106271_wp, 0.133081_wp, 0.135211_wp, &
                      0.133086_wp]
         case ("chelpg")
            charges = [ &
                      0.559048_wp, -0.247099_wp, 0.12256_wp, 0.143912_wp, -0.40774_wp, &
                      -0.008315_wp, -0.579405_wp, 0.072462_wp, 0.07239_wp, -0.000447_wp, &
                      0.000807_wp, -0.008273_wp, -0.007294_wp, 0.091338_wp, 0.104849_wp, &
                      0.091206_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.170554_wp, -0.032365_wp, -0.036063_wp, -0.031863_wp, -0.071411_wp, &
                      0.047874_wp, -0.313127_wp, 0.051527_wp, 0.051401_wp, 0.022065_wp, &
                      0.022401_wp, 0.023303_wp, 0.023351_wp, 0.02349_wp, 0.025338_wp, &
                      0.023529_wp]
         case ("resp")
            charges = [ &
                      0.555483_wp, -0.226456_wp, 0.112038_wp, 0.143014_wp, -0.395416_wp, &
                      -0.003659_wp, -0.589584_wp, 0.070631_wp, 0.070551_wp, -0.000525_wp, &
                      0.000801_wp, -0.007518_wp, -0.006536_wp, 0.088057_wp, 0.101174_wp, &
                      0.087945_wp]
         case ("mbis")
            charges = [ &
                      0.503702_wp, -0.41387_wp, -0.170887_wp, -0.16678_wp, -0.437496_wp, &
                      0.044658_wp, -0.54867_wp, 0.172633_wp, 0.172383_wp, 0.113893_wp, &
                      0.114027_wp, 0.10724_wp, 0.107222_wp, 0.1331_wp, 0.135721_wp, &
                      0.133128_wp]
         case ("chelpg")
            charges = [ &
                      0.565226_wp, -0.246685_wp, 0.12207_wp, 0.149873_wp, -0.412448_wp, &
                      -0.005557_wp, -0.591829_wp, 0.074765_wp, 0.074692_wp, -0.002153_wp, &
                      -0.000825_wp, -0.00863_wp, -0.007652_wp, 0.09179_wp, 0.105683_wp, &
                      0.091682_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (152)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.036595_wp, -0.034202_wp, -0.034202_wp, 0.019152_wp, 0.019153_wp, &
                      0.018247_wp, 0.018248_wp, 0.019153_wp, 0.019152_wp, -0.070304_wp, &
                      0.021475_wp, 0.021474_wp, 0.023304_wp, -0.070304_wp, 0.021475_wp, &
                      0.021474_wp, 0.023304_wp]
         case ("resp")
            charges = [ &
                      -0.056612_wp, 0.15672_wp, 0.161829_wp, -0.02585_wp, -0.025791_wp, &
                      0.01509_wp, 0.016305_wp, -0.027188_wp, -0.027626_wp, -0.289719_wp, &
                      0.062901_wp, 0.06351_wp, 0.070446_wp, -0.291558_wp, 0.063589_wp, &
                      0.062942_wp, 0.071013_wp]
         case ("mbis")
            charges = [ &
                      -0.213051_wp, -0.163978_wp, -0.163977_wp, 0.100155_wp, 0.100143_wp, &
                      0.106246_wp, 0.10624_wp, 0.100157_wp, 0.100142_wp, -0.420345_wp, &
                      0.127359_wp, 0.127348_wp, 0.1296_wp, -0.420346_wp, 0.127357_wp, &
                      0.127351_wp, 0.1296_wp]
         case ("chelpg")
            charges = [ &
                      -0.07034_wp, 0.17129_wp, 0.176535_wp, -0.028601_wp, -0.028541_wp, &
                      0.017467_wp, 0.018668_wp, -0.029973_wp, -0.03043_wp, -0.303698_wp, &
                      0.066158_wp, 0.066766_wp, 0.073322_wp, -0.305525_wp, 0.066839_wp, &
                      0.066178_wp, 0.073882_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.03696_wp, -0.034601_wp, -0.034601_wp, 0.019506_wp, 0.019508_wp, &
                      0.018622_wp, 0.018623_wp, 0.019507_wp, 0.019507_wp, -0.071523_wp, &
                      0.021721_wp, 0.021721_wp, 0.023527_wp, -0.071523_wp, 0.021721_wp, &
                      0.021721_wp, 0.023527_wp]
         case ("resp")
            charges = [ &
                      -0.056535_wp, 0.163906_wp, 0.16907_wp, -0.027432_wp, -0.027363_wp, &
                      0.014843_wp, 0.016073_wp, -0.028772_wp, -0.02922_wp, -0.297875_wp, &
                      0.064281_wp, 0.064892_wp, 0.072003_wp, -0.299767_wp, 0.064978_wp, &
                      0.064334_wp, 0.072584_wp]
         case ("mbis")
            charges = [ &
                      -0.214925_wp, -0.164309_wp, -0.164308_wp, 0.101027_wp, 0.101017_wp, &
                      0.107332_wp, 0.107326_wp, 0.101029_wp, 0.101016_wp, -0.426964_wp, &
                      0.129053_wp, 0.129044_wp, 0.131267_wp, -0.426964_wp, 0.129051_wp, &
                      0.129046_wp, 0.131267_wp]
         case ("chelpg")
            charges = [ &
                      -0.070368_wp, 0.178638_wp, 0.183933_wp, -0.030219_wp, -0.03015_wp, &
                      0.01723_wp, 0.018446_wp, -0.031592_wp, -0.032059_wp, -0.311973_wp, &
                      0.067563_wp, 0.068174_wp, 0.074903_wp, -0.31385_wp, 0.068254_wp, &
                      0.067595_wp, 0.075476_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.037384_wp, -0.035058_wp, -0.035058_wp, 0.019983_wp, 0.019986_wp, &
                      0.019109_wp, 0.01911_wp, 0.019984_wp, 0.019986_wp, -0.073061_wp, &
                      0.021995_wp, 0.021996_wp, 0.023743_wp, -0.073061_wp, 0.021995_wp, &
                      0.021996_wp, 0.023743_wp]
         case ("resp")
            charges = [ &
                      -0.05609_wp, 0.173741_wp, 0.178973_wp, -0.02962_wp, -0.029539_wp, &
                      0.01436_wp, 0.01561_wp, -0.030963_wp, -0.031422_wp, -0.308187_wp, &
                      0.065912_wp, 0.066526_wp, 0.07382_wp, -0.310146_wp, 0.066622_wp, &
                      0.065982_wp, 0.074418_wp]
         case ("mbis")
            charges = [ &
                      -0.217203_wp, -0.164626_wp, -0.164625_wp, 0.102149_wp, 0.10214_wp, &
                      0.108694_wp, 0.108688_wp, 0.10215_wp, 0.10214_wp, -0.435227_wp, &
                      0.131119_wp, 0.131113_wp, 0.133244_wp, -0.435228_wp, 0.131118_wp, &
                      0.131114_wp, 0.133244_wp]
         case ("chelpg")
            charges = [ &
                      -0.070013_wp, 0.188677_wp, 0.194035_wp, -0.032456_wp, -0.032375_wp, &
                      0.016747_wp, 0.017982_wp, -0.033831_wp, -0.034309_wp, -0.322431_wp, &
                      0.069225_wp, 0.069838_wp, 0.076748_wp, -0.324371_wp, 0.069927_wp, &
                      0.069271_wp, 0.077336_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (153)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068891_wp, -0.033585_wp, -0.034167_wp, -0.035421_wp, 0.233285_wp, &
                      -0.298093_wp, -0.18815_wp, 0.024296_wp, 0.024744_wp, 0.024259_wp, &
                      0.020259_wp, 0.020215_wp, 0.026195_wp, 0.026292_wp, 0.040617_wp, &
                      0.040291_wp, 0.177855_wp]
         case ("resp")
            charges = [ &
                      -0.3263_wp, 0.133868_wp, -0.007154_wp, -0.230988_wp, 0.726838_wp, &
                      -0.573169_wp, -0.609643_wp, 0.075956_wp, 0.083196_wp, 0.076814_wp, &
                      -0.010022_wp, -0.011057_wp, 0.036762_wp, 0.038433_wp, 0.085673_wp, &
                      0.087814_wp, 0.422977_wp]
         case ("mbis")
            charges = [ &
                      -0.421836_wp, -0.167217_wp, -0.17287_wp, -0.454703_wp, 0.844816_wp, &
                      -0.605428_wp, -0.657195_wp, 0.131139_wp, 0.131366_wp, 0.131066_wp, &
                      0.101489_wp, 0.101447_wp, 0.120962_wp, 0.121523_wp, 0.166943_wp, &
                      0.166934_wp, 0.461564_wp]
         case ("chelpg")
            charges = [ &
                      -0.344595_wp, 0.144491_wp, -0.005227_wp, -0.249549_wp, 0.736025_wp, &
                      -0.575156_wp, -0.61084_wp, 0.080021_wp, 0.087715_wp, 0.080873_wp, &
                      -0.011591_wp, -0.012616_wp, 0.037126_wp, 0.038802_wp, 0.089912_wp, &
                      0.092057_wp, 0.422552_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070563_wp, -0.032709_wp, -0.034968_wp, -0.034865_wp, 0.237532_wp, &
                      -0.321737_wp, -0.189025_wp, 0.024233_wp, 0.025007_wp, 0.024226_wp, &
                      0.02232_wp, 0.02227_wp, 0.024684_wp, 0.024613_wp, 0.045465_wp, &
                      0.045041_wp, 0.188475_wp]
         case ("resp")
            charges = [ &
                      -0.357033_wp, 0.159798_wp, 0.020919_wp, -0.262266_wp, 0.776101_wp, &
                      -0.623227_wp, -0.636418_wp, 0.081108_wp, 0.088919_wp, 0.082015_wp, &
                      -0.01272_wp, -0.013785_wp, 0.025495_wp, 0.02684_wp, 0.096196_wp, &
                      0.098277_wp, 0.449782_wp]
         case ("mbis")
            charges = [ &
                      -0.431982_wp, -0.165379_wp, -0.168567_wp, -0.473474_wp, 0.875855_wp, &
                      -0.644853_wp, -0.676683_wp, 0.132916_wp, 0.133737_wp, 0.132896_wp, &
                      0.104669_wp, 0.104616_wp, 0.118104_wp, 0.11836_wp, 0.177639_wp, &
                      0.177458_wp, 0.484687_wp]
         case ("chelpg")
            charges = [ &
                      -0.375664_wp, 0.169642_wp, 0.026341_wp, -0.283308_wp, 0.785849_wp, &
                      -0.625309_wp, -0.637419_wp, 0.085208_wp, 0.093632_wp, 0.086109_wp, &
                      -0.014235_wp, -0.015287_wp, 0.025062_wp, 0.02641_wp, 0.100842_wp, &
                      0.102922_wp, 0.449208_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071693_wp, -0.03186_wp, -0.03563_wp, -0.034296_wp, 0.240275_wp, &
                      -0.338734_wp, -0.188481_wp, 0.02388_wp, 0.025288_wp, 0.023905_wp, &
                      0.02417_wp, 0.024125_wp, 0.023001_wp, 0.02279_wp, 0.049213_wp, &
                      0.048735_wp, 0.195311_wp]
         case ("resp")
            charges = [ &
                      -0.379294_wp, 0.180769_wp, 0.038721_wp, -0.281775_wp, 0.808349_wp, &
                      -0.658621_wp, -0.651832_wp, 0.084328_wp, 0.093229_wp, 0.085286_wp, &
                      -0.014442_wp, -0.015504_wp, 0.01657_wp, 0.01763_wp, 0.103793_wp, &
                      0.105847_wp, 0.466945_wp]
         case ("mbis")
            charges = [ &
                      -0.438404_wp, -0.164221_wp, -0.164696_wp, -0.486328_wp, 0.896226_wp, &
                      -0.672803_wp, -0.686979_wp, 0.133517_wp, 0.135488_wp, 0.133554_wp, &
                      0.107559_wp, 0.107515_wp, 0.114903_wp, 0.114903_wp, 0.185601_wp, &
                      0.185323_wp, 0.498843_wp]
         case ("chelpg")
            charges = [ &
                      -0.398126_wp, 0.190159_wp, 0.046252_wp, -0.304273_wp, 0.81842_wp, &
                      -0.660755_wp, -0.652714_wp, 0.088448_wp, 0.098057_wp, 0.0894_wp, &
                      -0.01593_wp, -0.016978_wp, 0.015649_wp, 0.016712_wp, 0.108674_wp, &
                      0.110723_wp, 0.466282_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (154)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.032339_wp, -0.032995_wp, -0.03774_wp, 0.052123_wp, -0.148416_wp, &
                      0.229962_wp, -0.06926_wp, -0.303263_wp, 0.021116_wp, 0.021092_wp, &
                      0.021919_wp, 0.021946_wp, 0.02602_wp, 0.025867_wp, 0.032558_wp, &
                      0.032868_wp, 0.025064_wp, 0.022826_wp, 0.022835_wp, -0.066056_wp, &
                      0.045522_wp, 0.047677_wp, 0.040675_wp]
         case ("resp")
            charges = [ &
                      0.255686_wp, -0.071814_wp, -0.140534_wp, 0.281044_wp, -0.495267_wp, &
                      0.894584_wp, -0.365624_wp, -0.60234_wp, -0.038099_wp, -0.037376_wp, &
                      0.0183_wp, 0.019545_wp, 0.05396_wp, 0.053475_wp, 0.020175_wp, &
                      0.020047_wp, 0.084665_wp, 0.082499_wp, 0.082789_wp, -0.666288_wp, &
                      0.184283_wp, 0.184431_wp, 0.181856_wp]
         case ("mbis")
            charges = [ &
                      -0.161583_wp, -0.199354_wp, -0.294559_wp, 0.068779_wp, -0.413469_wp, &
                      0.842974_wp, -0.422741_wp, -0.60148_wp, 0.103038_wp, 0.102907_wp, &
                      0.109069_wp, 0.109069_wp, 0.130455_wp, 0.130106_wp, 0.095692_wp, &
                      0.095782_wp, 0.132413_wp, 0.129349_wp, 0.129339_wp, -0.66632_wp, &
                      0.197982_wp, 0.192598_wp, 0.189957_wp]
         case ("chelpg")
            charges = [ &
                      0.273047_wp, -0.080192_wp, -0.152604_wp, 0.30024_wp, -0.501884_wp, &
                      0.90283_wp, -0.385363_wp, -0.604511_wp, -0.040972_wp, -0.040282_wp, &
                      0.020255_wp, 0.021519_wp, 0.056618_wp, 0.056242_wp, 0.015885_wp, &
                      0.015781_wp, 0.088952_wp, 0.086918_wp, 0.087199_wp, -0.682414_wp, &
                      0.188124_wp, 0.188267_wp, 0.186343_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.032521_wp, -0.03194_wp, -0.037489_wp, 0.051441_wp, -0.148836_wp, &
                      0.234629_wp, -0.071358_wp, -0.335646_wp, 0.022026_wp, 0.02202_wp, &
                      0.024059_wp, 0.024065_wp, 0.027293_wp, 0.027214_wp, 0.034524_wp, &
                      0.034786_wp, 0.024968_wp, 0.023042_wp, 0.023047_wp, -0.065388_wp, &
                      0.048219_wp, 0.055252_wp, 0.046596_wp]
         case ("resp")
            charges = [ &
                      0.27602_wp, -0.062953_wp, -0.139545_wp, 0.313504_wp, -0.526944_wp, &
                      0.944782_wp, -0.388827_wp, -0.668289_wp, -0.041848_wp, -0.041103_wp, &
                      0.016845_wp, 0.018212_wp, 0.052483_wp, 0.052271_wp, 0.017897_wp, &
                      0.017375_wp, 0.088156_wp, 0.086391_wp, 0.086632_wp, -0.688864_wp, &
                      0.192074_wp, 0.199619_wp, 0.196111_wp]
         case ("mbis")
            charges = [ &
                      -0.16226_wp, -0.199299_wp, -0.296211_wp, 0.064406_wp, -0.421224_wp, &
                      0.878569_wp, -0.434108_wp, -0.653849_wp, 0.10485_wp, 0.104752_wp, &
                      0.112512_wp, 0.112477_wp, 0.132795_wp, 0.132579_wp, 0.100328_wp, &
                      0.100289_wp, 0.134531_wp, 0.131925_wp, 0.13191_wp, -0.691579_wp, &
                      0.205411_wp, 0.208311_wp, 0.202886_wp]
         case ("chelpg")
            charges = [ &
                      0.293427_wp, -0.070854_wp, -0.151919_wp, 0.332964_wp, -0.53357_wp, &
                      0.953002_wp, -0.408765_wp, -0.670452_wp, -0.04475_wp, -0.044037_wp, &
                      0.018672_wp, 0.02006_wp, 0.055168_wp, 0.055065_wp, 0.013541_wp, &
                      0.013043_wp, 0.092502_wp, 0.09085_wp, 0.09108_wp, -0.705044_wp, &
                      0.195934_wp, 0.20347_wp, 0.200614_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.03251_wp, -0.031554_wp, -0.037286_wp, 0.051231_wp, -0.148546_wp, &
                      0.235921_wp, -0.071971_wp, -0.346171_wp, 0.022326_wp, 0.022332_wp, &
                      0.024723_wp, 0.024715_wp, 0.027755_wp, 0.027717_wp, 0.034941_wp, &
                      0.035169_wp, 0.024842_wp, 0.023063_wp, 0.023065_wp, -0.065051_wp, &
                      0.048862_wp, 0.057732_wp, 0.048696_wp]
         case ("resp")
            charges = [ &
                      0.282595_wp, -0.060128_wp, -0.138635_wp, 0.323301_wp, -0.535401_wp, &
                      0.958906_wp, -0.395797_wp, -0.689374_wp, -0.043053_wp, -0.042299_wp, &
                      0.016344_wp, 0.017763_wp, 0.051978_wp, 0.051925_wp, 0.016944_wp, &
                      0.016207_wp, 0.089008_wp, 0.087488_wp, 0.087697_wp, -0.694481_wp, &
                      0.193866_wp, 0.204382_wp, 0.200763_wp]
         case ("mbis")
            charges = [ &
                      -0.162347_wp, -0.199196_wp, -0.296433_wp, 0.063256_wp, -0.422823_wp, &
                      0.88912_wp, -0.437314_wp, -0.67071_wp, 0.105402_wp, 0.105326_wp, &
                      0.113548_wp, 0.113489_wp, 0.133554_wp, 0.133412_wp, 0.101398_wp, &
                      0.101282_wp, 0.134967_wp, 0.132599_wp, 0.132577_wp, -0.698902_wp, &
                      0.207275_wp, 0.213294_wp, 0.207229_wp]
         case ("chelpg")
            charges = [ &
                      0.300011_wp, -0.067871_wp, -0.151106_wp, 0.342852_wp, -0.542037_wp, &
                      0.967122_wp, -0.415789_wp, -0.691537_wp, -0.045963_wp, -0.045242_wp, &
                      0.018129_wp, 0.019569_wp, 0.054671_wp, 0.054727_wp, 0.012566_wp, &
                      0.011852_wp, 0.093371_wp, 0.091956_wp, 0.092155_wp, -0.710683_wp, &
                      0.197732_wp, 0.20824_wp, 0.205274_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (155)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.00977_wp, -0.044231_wp, -0.233397_wp, 0.016393_wp, 0.016414_wp, &
                      -0.035776_wp, 0.024024_wp, 0.024016_wp, 0.090481_wp, 0.090473_wp, &
                      -0.033784_wp, 0.019554_wp, 0.019562_wp, -0.070018_wp, 0.019368_wp, &
                      0.019366_wp, 0.023868_wp, 0.02196_wp, 0.021959_wp]
         case ("resp")
            charges = [ &
                      0.511211_wp, -0.237702_wp, -0.960075_wp, 0.024523_wp, 0.023218_wp, &
                      0.043745_wp, -0.037265_wp, -0.039618_wp, 0.331863_wp, 0.331357_wp, &
                      0.302544_wp, -0.024943_wp, -0.023221_wp, -0.416788_wp, -0.050595_wp, &
                      -0.051183_wp, 0.091338_wp, 0.091406_wp, 0.090185_wp]
         case ("mbis")
            charges = [ &
                      0.067437_wp, -0.320778_wp, -0.831598_wp, 0.114665_wp, 0.114679_wp, &
                      -0.186647_wp, 0.08341_wp, 0.083397_wp, 0.331211_wp, 0.331177_wp, &
                      -0.160884_wp, 0.104201_wp, 0.104252_wp, -0.423673_wp, 0.100726_wp, &
                      0.100709_wp, 0.130872_wp, 0.128429_wp, 0.128418_wp]
         case ("chelpg")
            charges = [ &
                      0.528062_wp, -0.257239_wp, -0.965902_wp, 0.027823_wp, 0.026563_wp, &
                      0.052091_wp, -0.040457_wp, -0.042831_wp, 0.333123_wp, 0.332621_wp, &
                      0.314257_wp, -0.026744_wp, -0.025021_wp, -0.436845_wp, -0.052757_wp, &
                      -0.053312_wp, 0.096284_wp, 0.095742_wp, 0.094543_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.008089_wp, -0.042787_wp, -0.252865_wp, 0.019208_wp, 0.019224_wp, &
                      -0.035973_wp, 0.02368_wp, 0.023657_wp, 0.097492_wp, 0.097486_wp, &
                      -0.034291_wp, 0.020339_wp, 0.020339_wp, -0.072145_wp, 0.020059_wp, &
                      0.020057_wp, 0.023999_wp, 0.022218_wp, 0.022218_wp]
         case ("resp")
            charges = [ &
                      0.56298_wp, -0.208692_wp, -1.086268_wp, 0.022737_wp, 0.021391_wp, &
                      0.008478_wp, -0.047215_wp, -0.049754_wp, 0.374376_wp, 0.373848_wp, &
                      0.331461_wp, -0.019576_wp, -0.017841_wp, -0.425696_wp, -0.057083_wp, &
                      -0.057661_wp, 0.090738_wp, 0.09251_wp, 0.091266_wp]
         case ("mbis")
            charges = [ &
                      0.083755_wp, -0.318298_wp, -0.912734_wp, 0.118689_wp, 0.118687_wp, &
                      -0.193078_wp, 0.082176_wp, 0.082137_wp, 0.360694_wp, 0.360669_wp, &
                      -0.16068_wp, 0.106475_wp, 0.106511_wp, -0.434676_wp, 0.102165_wp, &
                      0.102151_wp, 0.133268_wp, 0.131051_wp, 0.131041_wp]
         case ("chelpg")
            charges = [ &
                      0.579417_wp, -0.226285_wp, -1.091894_wp, 0.025732_wp, 0.024429_wp, &
                      0.013463_wp, -0.050368_wp, -0.052929_wp, 0.375533_wp, 0.37501_wp, &
                      0.344653_wp, -0.020623_wp, -0.018886_wp, -0.446126_wp, -0.059452_wp, &
                      -0.059997_wp, 0.095653_wp, 0.096947_wp, 0.095725_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.007437_wp, -0.042029_wp, -0.260191_wp, 0.020527_wp, 0.020537_wp, &
                      -0.036018_wp, 0.023288_wp, 0.023258_wp, 0.100127_wp, 0.100122_wp, &
                      -0.034385_wp, 0.020524_wp, 0.02052_wp, -0.07289_wp, 0.020373_wp, &
                      0.020373_wp, 0.023969_wp, 0.022229_wp, 0.02223_wp]
         case ("resp")
            charges = [ &
                      0.581823_wp, -0.195328_wp, -1.135021_wp, 0.022208_wp, 0.020842_wp, &
                      -0.006598_wp, -0.051266_wp, -0.05388_wp, 0.390869_wp, 0.390332_wp, &
                      0.342609_wp, -0.017396_wp, -0.015655_wp, -0.428255_wp, -0.059479_wp, &
                      -0.060052_wp, 0.09018_wp, 0.092661_wp, 0.091406_wp]
         case ("mbis")
            charges = [ &
                      0.090338_wp, -0.317041_wp, -0.943133_wp, 0.120562_wp, 0.120551_wp, &
                      -0.195347_wp, 0.081207_wp, 0.081156_wp, 0.371642_wp, 0.37162_wp, &
                      -0.160489_wp, 0.107091_wp, 0.107119_wp, -0.438324_wp, 0.102748_wp, &
                      0.102736_wp, 0.133953_wp, 0.131811_wp, 0.131802_wp]
         case ("chelpg")
            charges = [ &
                      0.598047_wp, -0.211987_wp, -1.140576_wp, 0.025058_wp, 0.023734_wp, &
                      -0.003201_wp, -0.054395_wp, -0.057031_wp, 0.391987_wp, 0.391455_wp, &
                      0.356481_wp, -0.018086_wp, -0.016343_wp, -0.448826_wp, -0.061944_wp, &
                      -0.062485_wp, 0.09507_wp, 0.097137_wp, 0.095905_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (156)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.081961_wp, 0.082082_wp, 0.082068_wp, 0.081962_wp, 0.082082_wp, &
                      0.082069_wp, -0.082035_wp, -0.082061_wp, -0.082032_wp, -0.082061_wp, &
                      -0.08203_wp, -0.082038_wp]
         case ("resp")
            charges = [ &
                      0.101572_wp, 0.09906_wp, 0.098842_wp, 0.103101_wp, 0.096439_wp, &
                      0.104013_wp, -0.100987_wp, -0.099838_wp, -0.099894_wp, -0.099671_wp, &
                      -0.101366_wp, -0.101271_wp]
         case ("mbis")
            charges = [ &
                      0.131697_wp, 0.131901_wp, 0.131879_wp, 0.131701_wp, 0.131897_wp, &
                      0.131883_wp, -0.131812_wp, -0.131859_wp, -0.131825_wp, -0.131859_wp, &
                      -0.131822_wp, -0.131815_wp]
         case ("chelpg")
            charges = [ &
                      0.101462_wp, 0.099313_wp, 0.098709_wp, 0.103439_wp, 0.096242_wp, &
                      0.104336_wp, -0.10111_wp, -0.099936_wp, -0.099942_wp, -0.099701_wp, &
                      -0.101485_wp, -0.101327_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.085332_wp, 0.0855_wp, 0.085483_wp, 0.085333_wp, 0.0855_wp, &
                      0.085484_wp, -0.085457_wp, -0.085452_wp, -0.085423_wp, -0.085453_wp, &
                      -0.08542_wp, -0.08546_wp]
         case ("resp")
            charges = [ &
                      0.108541_wp, 0.106464_wp, 0.106227_wp, 0.110028_wp, 0.10382_wp, &
                      0.11147_wp, -0.108196_wp, -0.107099_wp, -0.107151_wp, -0.106932_wp, &
                      -0.10865_wp, -0.108521_wp]
         case ("mbis")
            charges = [ &
                      0.137328_wp, 0.137683_wp, 0.13765_wp, 0.137332_wp, 0.137679_wp, &
                      0.137654_wp, -0.13755_wp, -0.137582_wp, -0.137546_wp, -0.137582_wp, &
                      -0.137544_wp, -0.137554_wp]
         case ("chelpg")
            charges = [ &
                      0.108447_wp, 0.106705_wp, 0.10612_wp, 0.110338_wp, 0.103655_wp, &
                      0.111773_wp, -0.108316_wp, -0.107199_wp, -0.107205_wp, -0.10697_wp, &
                      -0.108768_wp, -0.108579_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.089072_wp, 0.08929_wp, 0.089269_wp, 0.089073_wp, 0.08929_wp, &
                      0.089271_wp, -0.089254_wp, -0.089212_wp, -0.089182_wp, -0.089212_wp, &
                      -0.089179_wp, -0.089258_wp]
         case ("resp")
            charges = [ &
                      0.116311_wp, 0.114702_wp, 0.114448_wp, 0.117747_wp, 0.112037_wp, &
                      0.119763_wp, -0.116229_wp, -0.115181_wp, -0.115228_wp, -0.115015_wp, &
                      -0.116757_wp, -0.116598_wp]
         case ("mbis")
            charges = [ &
                      0.143584_wp, 0.144098_wp, 0.144053_wp, 0.143588_wp, 0.144094_wp, &
                      0.144057_wp, -0.143924_wp, -0.143934_wp, -0.143896_wp, -0.143934_wp, &
                      -0.143893_wp, -0.143927_wp]
         case ("chelpg")
            charges = [ &
                      0.116234_wp, 0.114929_wp, 0.114367_wp, 0.118029_wp, 0.111904_wp, &
                      0.120047_wp, -0.116346_wp, -0.115282_wp, -0.115289_wp, -0.115061_wp, &
                      -0.116874_wp, -0.11666_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (157)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.257477_wp, 0.15713_wp, 0.031413_wp, -0.034251_wp, 0.039706_wp, &
                      0.039689_wp, -0.065876_wp, 0.027388_wp, 0.031139_wp, 0.031133_wp]
         case ("resp")
            charges = [ &
                      -0.448536_wp, 0.384273_wp, -0.016764_wp, 0.119592_wp, 0.000265_wp, &
                      -0.000469_wp, -0.273749_wp, 0.068594_wp, 0.083022_wp, 0.083772_wp]
         case ("mbis")
            charges = [ &
                      -0.456158_wp, 0.448395_wp, 0.021314_wp, -0.330224_wp, 0.142505_wp, &
                      0.142473_wp, -0.384073_wp, 0.128994_wp, 0.143398_wp, 0.143372_wp]
         case ("chelpg")
            charges = [ &
                      -0.44865_wp, 0.383975_wp, -0.01722_wp, 0.126836_wp, -0.001015_wp, &
                      -0.001748_wp, -0.287348_wp, 0.072004_wp, 0.086211_wp, 0.086956_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.309285_wp, 0.169659_wp, 0.047227_wp, -0.029708_wp, 0.052046_wp, &
                      0.052031_wp, -0.069088_wp, 0.032983_wp, 0.027062_wp, 0.027068_wp]
         case ("resp")
            charges = [ &
                      -0.540151_wp, 0.426152_wp, 0.003586_wp, 0.137903_wp, 0.013536_wp, &
                      0.012866_wp, -0.293551_wp, 0.083597_wp, 0.077666_wp, 0.078396_wp]
         case ("mbis")
            charges = [ &
                      -0.537213_wp, 0.489319_wp, 0.043729_wp, -0.349657_wp, 0.164387_wp, &
                      0.164358_wp, -0.392642_wp, 0.141164_wp, 0.138277_wp, 0.138272_wp]
         case ("chelpg")
            charges = [ &
                      -0.540233_wp, 0.425774_wp, 0.00316_wp, 0.145146_wp, 0.012273_wp, &
                      0.011604_wp, -0.307266_wp, 0.087037_wp, 0.08089_wp, 0.081616_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.313061_wp, 0.170525_wp, 0.048408_wp, -0.029346_wp, 0.053017_wp, &
                      0.053002_wp, -0.069312_wp, 0.033463_wp, 0.026646_wp, 0.026653_wp]
         case ("resp")
            charges = [ &
                      -0.546813_wp, 0.428991_wp, 0.005189_wp, 0.139202_wp, 0.014661_wp, &
                      0.013996_wp, -0.294839_wp, 0.084871_wp, 0.077007_wp, 0.077736_wp]
         case ("mbis")
            charges = [ &
                      -0.543115_wp, 0.492165_wp, 0.045426_wp, -0.35108_wp, 0.166091_wp, &
                      0.166061_wp, -0.393015_wp, 0.142154_wp, 0.137656_wp, 0.137653_wp]
         case ("chelpg")
            charges = [ &
                      -0.546894_wp, 0.428609_wp, 0.004764_wp, 0.146448_wp, 0.013398_wp, &
                      0.012734_wp, -0.308558_wp, 0.088312_wp, 0.080231_wp, 0.080955_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (158)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.23328_wp, -0.299451_wp, -0.190663_wp, 0.178822_wp, -0.031349_wp, &
                      -0.065122_wp, 0.042388_wp, 0.042053_wp, 0.028176_wp, 0.030982_wp, &
                      0.030881_wp]
         case ("resp")
            charges = [ &
                      0.62501_wp, -0.538022_wp, -0.625251_wp, 0.436751_wp, 0.091232_wp, &
                      -0.399527_wp, 0.03972_wp, 0.039895_wp, 0.098436_wp, 0.115863_wp, &
                      0.115892_wp]
         case ("mbis")
            charges = [ &
                      0.824126_wp, -0.596079_wp, -0.65508_wp, 0.460642_wp, -0.391511_wp, &
                      -0.380963_wp, 0.161191_wp, 0.161237_wp, 0.130529_wp, 0.143329_wp, &
                      0.142578_wp]
         case ("chelpg")
            charges = [ &
                      0.625588_wp, -0.538577_wp, -0.62707_wp, 0.437599_wp, 0.098573_wp, &
                      -0.415737_wp, 0.038684_wp, 0.038886_wp, 0.102429_wp, 0.119809_wp, &
                      0.119817_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.238188_wp, -0.327109_wp, -0.191738_wp, 0.190991_wp, -0.03027_wp, &
                      -0.066704_wp, 0.04829_wp, 0.047827_wp, 0.032521_wp, 0.02891_wp, &
                      0.029093_wp]
         case ("resp")
            charges = [ &
                      0.668326_wp, -0.592835_wp, -0.656502_wp, 0.469018_wp, 0.087494_wp, &
                      -0.391945_wp, 0.047329_wp, 0.047277_wp, 0.106746_wp, 0.107233_wp, &
                      0.107858_wp]
         case ("mbis")
            charges = [ &
                      0.858248_wp, -0.641339_wp, -0.677023_wp, 0.486916_wp, -0.408733_wp, &
                      -0.384518_wp, 0.173244_wp, 0.17307_wp, 0.139266_wp, 0.140558_wp, &
                      0.14031_wp]
         case ("chelpg")
            charges = [ &
                      0.669009_wp, -0.593401_wp, -0.658322_wp, 0.469853_wp, 0.094498_wp, &
                      -0.407944_wp, 0.046373_wp, 0.046347_wp, 0.1107_wp, 0.111142_wp, &
                      0.111746_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.240363_wp, -0.340929_wp, -0.191316_wp, 0.196411_wp, -0.029591_wp, &
                      -0.06744_wp, 0.051414_wp, 0.050898_wp, 0.035068_wp, 0.027386_wp, &
                      0.027733_wp]
         case ("resp")
            charges = [ &
                      0.687791_wp, -0.619935_wp, -0.668935_wp, 0.483394_wp, 0.086255_wp, &
                      -0.387783_wp, 0.05142_wp, 0.051262_wp, 0.111786_wp, 0.101884_wp, &
                      0.102861_wp]
         case ("mbis")
            charges = [ &
                      0.873712_wp, -0.663641_wp, -0.684994_wp, 0.498057_wp, -0.416928_wp, &
                      -0.38539_wp, 0.17944_wp, 0.179173_wp, 0.144204_wp, 0.138162_wp, &
                      0.138206_wp]
         case ("chelpg")
            charges = [ &
                      0.688527_wp, -0.620509_wp, -0.670761_wp, 0.484227_wp, 0.093093_wp, &
                      -0.403671_wp, 0.050503_wp, 0.050373_wp, 0.115721_wp, 0.10577_wp, &
                      0.106726_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (159)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.248959_wp, 0.081545_wp, 0.001456_wp, 0.054368_wp, 0.054368_wp, &
                      -0.053462_wp, 0.035529_wp, 0.037579_wp, 0.037579_wp]
         case ("resp")
            charges = [ &
                      -0.443664_wp, 0.246137_wp, 0.094454_wp, 0.049341_wp, 0.048171_wp, &
                      -0.241959_wp, 0.06825_wp, 0.089757_wp, 0.089513_wp]
         case ("mbis")
            charges = [ &
                      -0.462137_wp, 0.362504_wp, -0.266212_wp, 0.161121_wp, 0.161122_wp, &
                      -0.397924_wp, 0.142023_wp, 0.14975_wp, 0.149755_wp]
         case ("chelpg")
            charges = [ &
                      -0.443583_wp, 0.245032_wp, 0.100752_wp, 0.048218_wp, 0.047045_wp, &
                      -0.253385_wp, 0.071026_wp, 0.092573_wp, 0.092322_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.320168_wp, 0.103556_wp, 0.008474_wp, 0.068782_wp, 0.068782_wp, &
                      -0.050858_wp, 0.043752_wp, 0.038842_wp, 0.03884_wp]
         case ("resp")
            charges = [ &
                      -0.560128_wp, 0.314894_wp, 0.091841_wp, 0.068585_wp, 0.067441_wp, &
                      -0.23342_wp, 0.081423_wp, 0.084794_wp, 0.084571_wp]
         case ("mbis")
            charges = [ &
                      -0.586805_wp, 0.45495_wp, -0.302223_wp, 0.187828_wp, 0.187828_wp, &
                      -0.399534_wp, 0.156172_wp, 0.150893_wp, 0.150892_wp]
         case ("chelpg")
            charges = [ &
                      -0.560103_wp, 0.314014_wp, 0.097506_wp, 0.067616_wp, 0.066471_wp, &
                      -0.244572_wp, 0.084169_wp, 0.087564_wp, 0.087337_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.323393_wp, 0.10454_wp, 0.008809_wp, 0.069447_wp, 0.069447_wp, &
                      -0.050722_wp, 0.044165_wp, 0.038854_wp, 0.038852_wp]
         case ("resp")
            charges = [ &
                      -0.565455_wp, 0.318059_wp, 0.091662_wp, 0.069506_wp, 0.068364_wp, &
                      -0.233008_wp, 0.082138_wp, 0.084478_wp, 0.084256_wp]
         case ("mbis")
            charges = [ &
                      -0.592436_wp, 0.459099_wp, -0.303841_wp, 0.189049_wp, 0.189049_wp, &
                      -0.39947_wp, 0.156865_wp, 0.150843_wp, 0.150842_wp]
         case ("chelpg")
            charges = [ &
                      -0.565436_wp, 0.317191_wp, 0.097294_wp, 0.068545_wp, 0.067402_wp, &
                      -0.244142_wp, 0.084881_wp, 0.087245_wp, 0.087019_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (160)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.230072_wp, -0.303219_wp, -0.148325_wp, -0.066013_wp, 0.04569_wp, &
                      0.040568_wp, 0.04783_wp, 0.051979_wp, -0.035105_wp, 0.033004_wp, &
                      0.032737_wp, -0.065242_wp, 0.027173_wp, 0.027307_wp, 0.029245_wp, &
                      0.026181_wp, 0.026114_wp]
         case ("resp")
            charges = [ &
                      0.922526_wp, -0.60631_wp, -0.528615_wp, -0.675334_wp, 0.18336_wp, &
                      0.186623_wp, 0.186777_wp, 0.210395_wp, 0.099794_wp, 0.037455_wp, &
                      0.036775_wp, -0.407055_wp, 0.023305_wp, 0.02445_wp, 0.112347_wp, &
                      0.096949_wp, 0.096557_wp]
         case ("mbis")
            charges = [ &
                      0.843825_wp, -0.60165_wp, -0.415695_wp, -0.666425_wp, 0.198082_wp, &
                      0.189966_wp, 0.192795_wp, 0.065042_wp, -0.240072_wp, 0.096935_wp, &
                      0.096271_wp, -0.408551_wp, 0.124569_wp, 0.12477_wp, 0.137401_wp, &
                      0.131426_wp, 0.131309_wp]
         case ("chelpg")
            charges = [ &
                      0.932364_wp, -0.608275_wp, -0.534098_wp, -0.693545_wp, 0.187539_wp, &
                      0.191553_wp, 0.191035_wp, 0.216691_wp, 0.107835_wp, 0.035373_wp, &
                      0.034669_wp, -0.424117_wp, 0.021977_wp, 0.023019_wp, 0.116543_wp, &
                      0.100923_wp, 0.100513_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.234929_wp, -0.337106_wp, -0.148653_wp, -0.065305_wp, 0.048516_wp, &
                      0.046697_wp, 0.05568_wp, 0.051307_wp, -0.034621_wp, 0.034938_wp, &
                      0.034917_wp, -0.06529_wp, 0.028619_wp, 0.028702_wp, 0.030749_wp, &
                      0.027974_wp, 0.027945_wp]
         case ("resp")
            charges = [ &
                      0.97696_wp, -0.675456_wp, -0.563042_wp, -0.700117_wp, 0.191547_wp, &
                      0.201932_wp, 0.202794_wp, 0.235259_wp, 0.121551_wp, 0.036679_wp, &
                      0.036547_wp, -0.422066_wp, 0.019524_wp, 0.020439_wp, 0.116836_wp, &
                      0.100513_wp, 0.100101_wp]
         case ("mbis")
            charges = [ &
                      0.880924_wp, -0.656417_wp, -0.423655_wp, -0.69266_wp, 0.205842_wp, &
                      0.203361_wp, 0.209085_wp, 0.060153_wp, -0.238654_wp, 0.101546_wp, &
                      0.101336_wp, -0.41609_wp, 0.126768_wp, 0.126878_wp, 0.140849_wp, &
                      0.135391_wp, 0.135338_wp]
         case ("chelpg")
            charges = [ &
                      0.986872_wp, -0.677416_wp, -0.568557_wp, -0.71847_wp, 0.195755_wp, &
                      0.206899_wp, 0.207083_wp, 0.241459_wp, 0.13019_wp, 0.034583_wp, &
                      0.034426_wp, -0.439532_wp, 0.018047_wp, 0.018859_wp, 0.121108_wp, &
                      0.104562_wp, 0.104132_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.236031_wp, -0.346178_wp, -0.148371_wp, -0.065008_wp, 0.049075_wp, &
                      0.048497_wp, 0.057797_wp, 0.05113_wp, -0.034404_wp, 0.035237_wp, &
                      0.035306_wp, -0.065239_wp, 0.029062_wp, 0.029113_wp, 0.031163_wp, &
                      0.028398_wp, 0.028389_wp]
         case ("resp")
            charges = [ &
                      0.989658_wp, -0.693681_wp, -0.570638_wp, -0.705353_wp, 0.19313_wp, &
                      0.206077_wp, 0.206951_wp, 0.241459_wp, 0.127387_wp, 0.036145_wp, &
                      0.036252_wp, -0.425547_wp, 0.01862_wp, 0.019408_wp, 0.117966_wp, &
                      0.101292_wp, 0.100873_wp]
         case ("mbis")
            charges = [ &
                      0.889981_wp, -0.670945_wp, -0.424987_wp, -0.698925_wp, 0.20745_wp, &
                      0.207084_wp, 0.21334_wp, 0.059087_wp, -0.238134_wp, 0.102356_wp, &
                      0.102312_wp, -0.417753_wp, 0.127399_wp, 0.127453_wp, 0.141724_wp, &
                      0.136285_wp, 0.13627_wp]
         case ("chelpg")
            charges = [ &
                      0.999589_wp, -0.695652_wp, -0.576159_wp, -0.723742_wp, 0.197346_wp, &
                      0.211052_wp, 0.21125_wp, 0.247645_wp, 0.136169_wp, 0.034046_wp, &
                      0.034128_wp, -0.443098_wp, 0.017106_wp, 0.01779_wp, 0.122253_wp, &
                      0.105356_wp, 0.10492_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (161)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.009713_wp, -0.042125_wp, -0.23361_wp, 0.017265_wp, 0.017271_wp, &
                      -0.06962_wp, 0.024251_wp, 0.024237_wp, 0.09069_wp, 0.090692_wp, &
                      0.024977_wp, 0.023132_wp, 0.023133_wp]
         case ("resp")
            charges = [ &
                      0.352022_wp, 0.047449_wp, -0.914306_wp, -0.009388_wp, -0.010723_wp, &
                      -0.299002_wp, -0.004329_wp, -0.005329_wp, 0.31837_wp, 0.317607_wp, &
                      0.092227_wp, 0.057248_wp, 0.058154_wp]
         case ("mbis")
            charges = [ &
                      0.061678_wp, -0.261481_wp, -0.83617_wp, 0.108298_wp, 0.108297_wp, &
                      -0.3983_wp, 0.084726_wp, 0.084717_wp, 0.332428_wp, 0.332427_wp, &
                      0.130376_wp, 0.126502_wp, 0.126503_wp]
         case ("chelpg")
            charges = [ &
                      0.356337_wp, 0.051932_wp, -0.915725_wp, -0.01036_wp, -0.011695_wp, &
                      -0.310316_wp, -0.005569_wp, -0.006579_wp, 0.318395_wp, 0.317625_wp, &
                      0.095275_wp, 0.059883_wp, 0.060797_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.008016_wp, -0.04048_wp, -0.253927_wp, 0.020319_wp, 0.020327_wp, &
                      -0.071239_wp, 0.023959_wp, 0.023932_wp, 0.09806_wp, 0.098062_wp, &
                      0.025593_wp, 0.02369_wp, 0.02369_wp]
         case ("resp")
            charges = [ &
                      0.411319_wp, 0.075886_wp, -1.049835_wp, -0.009889_wp, -0.011311_wp, &
                      -0.351407_wp, -0.015277_wp, -0.016417_wp, 0.363924_wp, 0.363074_wp, &
                      0.103304_wp, 0.067832_wp, 0.068795_wp]
         case ("mbis")
            charges = [ &
                      0.078901_wp, -0.256948_wp, -0.921501_wp, 0.112506_wp, 0.112508_wp, &
                      -0.413248_wp, 0.083546_wp, 0.083515_wp, 0.363471_wp, 0.363472_wp, &
                      0.133968_wp, 0.129907_wp, 0.129905_wp]
         case ("chelpg")
            charges = [ &
                      0.415482_wp, 0.081674_wp, -1.051128_wp, -0.011173_wp, -0.012597_wp, &
                      -0.363867_wp, -0.016565_wp, -0.017716_wp, 0.36389_wp, 0.363032_wp, &
                      0.106597_wp, 0.070699_wp, 0.071672_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.007485_wp, -0.039818_wp, -0.260049_wp, 0.021455_wp, 0.021463_wp, &
                      -0.071655_wp, 0.023647_wp, 0.023616_wp, 0.100285_wp, 0.100287_wp, &
                      0.025795_wp, 0.023747_wp, 0.023746_wp]
         case ("resp")
            charges = [ &
                      0.429271_wp, 0.085667_wp, -1.091958_wp, -0.009813_wp, -0.01126_wp, &
                      -0.368125_wp, -0.018967_wp, -0.020153_wp, 0.378139_wp, 0.377261_wp, &
                      0.106869_wp, 0.071045_wp, 0.072025_wp]
         case ("mbis")
            charges = [ &
                      0.084477_wp, -0.255446_wp, -0.947111_wp, 0.114089_wp, 0.114091_wp, &
                      -0.417339_wp, 0.082763_wp, 0.082724_wp, 0.372714_wp, 0.372716_wp, &
                      0.135008_wp, 0.130659_wp, 0.130656_wp]
         case ("chelpg")
            charges = [ &
                      0.433411_wp, 0.091802_wp, -1.093221_wp, -0.011182_wp, -0.012631_wp, &
                      -0.380889_wp, -0.020274_wp, -0.02147_wp, 0.378088_wp, 0.377203_wp, &
                      0.110226_wp, 0.073974_wp, 0.074963_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (162)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.19165_wp, 0.036796_wp, -0.039137_wp, -0.013947_wp, -0.039138_wp, &
                      0.036794_wp, 0.041397_wp, 0.041369_wp, 0.044751_wp, 0.04137_wp, &
                      0.041396_wp]
         case ("resp")
            charges = [ &
                      -0.64342_wp, 0.456533_wp, -0.511236_wp, 0.252094_wp, -0.513896_wp, &
                      0.458268_wp, 0.029602_wp, 0.183022_wp, 0.075696_wp, 0.184311_wp, &
                      0.029026_wp]
         case ("mbis")
            charges = [ &
                      -0.435643_wp, 0.173433_wp, -0.263683_wp, -0.030226_wp, -0.263686_wp, &
                      0.173426_wp, 0.099163_wp, 0.157886_wp, 0.132282_wp, 0.157886_wp, &
                      0.099163_wp]
         case ("chelpg")
            charges = [ &
                      -0.655405_wp, 0.475525_wp, -0.537999_wp, 0.27952_wp, -0.540509_wp, &
                      0.477221_wp, 0.026496_wp, 0.188813_wp, 0.070382_wp, 0.190038_wp, &
                      0.025918_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.229593_wp, 0.028849_wp, -0.037416_wp, -0.00921_wp, -0.037416_wp, &
                      0.028847_wp, 0.046804_wp, 0.052654_wp, 0.057024_wp, 0.052655_wp, &
                      0.046803_wp]
         case ("resp")
            charges = [ &
                      -0.75466_wp, 0.49263_wp, -0.541397_wp, 0.272571_wp, -0.544245_wp, &
                      0.494601_wp, 0.035882_wp, 0.207313_wp, 0.093435_wp, 0.208648_wp, &
                      0.035222_wp]
         case ("mbis")
            charges = [ &
                      -0.509344_wp, 0.173314_wp, -0.270906_wp, -0.030609_wp, -0.270908_wp, &
                      0.173306_wp, 0.112656_wp, 0.17837_wp, 0.153097_wp, 0.17837_wp, &
                      0.112656_wp]
         case ("chelpg")
            charges = [ &
                      -0.766654_wp, 0.511633_wp, -0.568181_wp, 0.299995_wp, -0.570879_wp, &
                      0.513565_wp, 0.032783_wp, 0.213108_wp, 0.088131_wp, 0.214379_wp, &
                      0.03212_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.233801_wp, 0.027865_wp, -0.03703_wp, -0.008407_wp, -0.03703_wp, &
                      0.027863_wp, 0.047178_wp, 0.053884_wp, 0.058417_wp, 0.053885_wp, &
                      0.047177_wp]
         case ("resp")
            charges = [ &
                      -0.766773_wp, 0.496175_wp, -0.544039_wp, 0.275291_wp, -0.546908_wp, &
                      0.498171_wp, 0.036278_wp, 0.209793_wp, 0.095271_wp, 0.211133_wp, &
                      0.035608_wp]
         case ("mbis")
            charges = [ &
                      -0.51752_wp, 0.173159_wp, -0.271208_wp, -0.030157_wp, -0.27121_wp, &
                      0.17315_wp, 0.113787_wp, 0.180463_wp, 0.155286_wp, 0.180463_wp, &
                      0.113787_wp]
         case ("chelpg")
            charges = [ &
                      -0.778774_wp, 0.515181_wp, -0.570828_wp, 0.302724_wp, -0.573547_wp, &
                      0.517139_wp, 0.033178_wp, 0.215589_wp, 0.089966_wp, 0.216865_wp, &
                      0.032506_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (163)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.050315_wp, 0.050316_wp, -0.025152_wp, -0.025155_wp, -0.025161_wp, &
                      -0.025164_wp]
         case ("resp")
            charges = [ &
                      0.006859_wp, 0.006606_wp, -0.003584_wp, -0.003259_wp, -0.003542_wp, &
                      -0.00308_wp]
         case ("mbis")
            charges = [ &
                      0.082453_wp, 0.082456_wp, -0.04122_wp, -0.041224_wp, -0.04123_wp, &
                      -0.041234_wp]
         case ("chelpg")
            charges = [ &
                      0.006885_wp, 0.006621_wp, -0.003596_wp, -0.003271_wp, -0.003551_wp, &
                      -0.003089_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.05125_wp, 0.05125_wp, -0.025619_wp, -0.025621_wp, -0.025628_wp, &
                      -0.025631_wp]
         case ("resp")
            charges = [ &
                      0.007539_wp, 0.007299_wp, -0.003928_wp, -0.003599_wp, -0.003893_wp, &
                      -0.003419_wp]
         case ("mbis")
            charges = [ &
                      0.084837_wp, 0.084839_wp, -0.042411_wp, -0.042415_wp, -0.042423_wp, &
                      -0.042427_wp]
         case ("chelpg")
            charges = [ &
                      0.007568_wp, 0.007318_wp, -0.003941_wp, -0.003612_wp, -0.003904_wp, &
                      -0.00343_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.052021_wp, 0.052021_wp, -0.026003_wp, -0.026006_wp, -0.026014_wp, &
                      -0.026017_wp]
         case ("resp")
            charges = [ &
                      0.008153_wp, 0.007924_wp, -0.004238_wp, -0.003905_wp, -0.004209_wp, &
                      -0.003726_wp]
         case ("mbis")
            charges = [ &
                      0.086805_wp, 0.086805_wp, -0.043394_wp, -0.043398_wp, -0.043407_wp, &
                      -0.043411_wp]
         case ("chelpg")
            charges = [ &
                      0.008184_wp, 0.007945_wp, -0.004252_wp, -0.003919_wp, -0.004221_wp, &
                      -0.003737_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (164)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.213216_wp, 0.040101_wp, -0.040423_wp, -0.040415_wp, 0.040112_wp, &
                      0.020896_wp, 0.030104_wp, 0.030104_wp, 0.020924_wp, 0.029425_wp, &
                      0.029391_wp, 0.0265_wp, 0.026489_wp]
         case ("resp")
            charges = [ &
                      -0.420018_wp, 0.146121_wp, -0.040861_wp, -0.050684_wp, 0.141607_wp, &
                      0.029074_wp, 0.007591_wp, 0.010881_wp, 0.031386_wp, 0.023931_wp, &
                      0.025465_wp, 0.046294_wp, 0.049213_wp]
         case ("mbis")
            charges = [ &
                      -0.396593_wp, 0.092002_wp, -0.286143_wp, -0.286107_wp, 0.092046_wp, &
                      0.062543_wp, 0.128156_wp, 0.12815_wp, 0.062589_wp, 0.076902_wp, &
                      0.076816_wp, 0.124837_wp, 0.124796_wp]
         case ("chelpg")
            charges = [ &
                      -0.422454_wp, 0.152082_wp, -0.045092_wp, -0.055374_wp, 0.147609_wp, &
                      0.027728_wp, 0.008884_wp, 0.012285_wp, 0.030057_wp, 0.022475_wp, &
                      0.024027_wp, 0.047373_wp, 0.050401_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.238229_wp, 0.037598_wp, -0.038642_wp, -0.038631_wp, 0.03761_wp, &
                      0.023715_wp, 0.035599_wp, 0.035596_wp, 0.023735_wp, 0.031128_wp, &
                      0.031099_wp, 0.029708_wp, 0.029708_wp]
         case ("resp")
            charges = [ &
                      -0.49367_wp, 0.168261_wp, -0.035076_wp, -0.044307_wp, 0.163416_wp, &
                      0.031441_wp, 0.012336_wp, 0.015492_wp, 0.033772_wp, 0.02349_wp, &
                      0.025148_wp, 0.048451_wp, 0.051247_wp]
         case ("mbis")
            charges = [ &
                      -0.438974_wp, 0.089202_wp, -0.287632_wp, -0.287598_wp, 0.089255_wp, &
                      0.069219_wp, 0.137191_wp, 0.137179_wp, 0.069248_wp, 0.081679_wp, &
                      0.0816_wp, 0.129823_wp, 0.1298_wp]
         case ("chelpg")
            charges = [ &
                      -0.496149_wp, 0.174309_wp, -0.03929_wp, -0.048978_wp, 0.169515_wp, &
                      0.030067_wp, 0.013619_wp, 0.016884_wp, 0.032411_wp, 0.022008_wp, &
                      0.023681_wp, 0.04951_wp, 0.052414_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.242983_wp, 0.037089_wp, -0.038197_wp, -0.038185_wp, 0.037102_wp, &
                      0.02419_wp, 0.036774_wp, 0.036771_wp, 0.024207_wp, 0.031307_wp, &
                      0.031279_wp, 0.030319_wp, 0.030321_wp]
         case ("resp")
            charges = [ &
                      -0.507827_wp, 0.172327_wp, -0.033561_wp, -0.042655_wp, 0.167407_wp, &
                      0.031803_wp, 0.013434_wp, 0.016557_wp, 0.034138_wp, 0.02323_wp, &
                      0.024918_wp, 0.048731_wp, 0.051498_wp]
         case ("mbis")
            charges = [ &
                      -0.447123_wp, 0.088761_wp, -0.287689_wp, -0.287656_wp, 0.088816_wp, &
                      0.070381_wp, 0.139073_wp, 0.13906_wp, 0.070406_wp, 0.08233_wp, &
                      0.082253_wp, 0.1307_wp, 0.130682_wp]
         case ("chelpg")
            charges = [ &
                      -0.510312_wp, 0.178382_wp, -0.037755_wp, -0.047307_wp, 0.173514_wp, &
                      0.030425_wp, 0.014709_wp, 0.017942_wp, 0.032773_wp, 0.021746_wp, &
                      0.023447_wp, 0.04978_wp, 0.052656_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (165)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.341311_wp, 0.440718_wp, -0.341291_wp, -0.028905_wp, -0.029975_wp, &
                      -0.029953_wp, -0.028926_wp, 0.055514_wp, 0.049613_wp, 0.03589_wp, &
                      0.038879_wp, 0.038869_wp, 0.035804_wp, 0.055413_wp, 0.049648_wp]
         case ("resp")
            charges = [ &
                      -0.526155_wp, 0.899347_wp, -0.525925_wp, -0.24262_wp, -0.078925_wp, &
                      -0.068812_wp, -0.233375_wp, 0.130162_wp, 0.126666_wp, 0.070504_wp, &
                      0.067671_wp, 0.063544_wp, 0.068083_wp, 0.126104_wp, 0.123731_wp]
         case ("mbis")
            charges = [ &
                      -0.618327_wp, 1.185761_wp, -0.618342_wp, -0.462728_wp, -0.183825_wp, &
                      -0.183378_wp, -0.462566_wp, 0.206529_wp, 0.199375_wp, 0.134368_wp, &
                      0.131723_wp, 0.131629_wp, 0.133993_wp, 0.206349_wp, 0.199425_wp]
         case ("chelpg")
            charges = [ &
                      -0.527709_wp, 0.906668_wp, -0.527462_wp, -0.250184_wp, -0.080937_wp, &
                      -0.070104_wp, -0.241085_wp, 0.132358_wp, 0.129027_wp, 0.072071_wp, &
                      0.068879_wp, 0.064578_wp, 0.069491_wp, 0.128309_wp, 0.126099_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.387456_wp, 0.4388_wp, -0.387431_wp, -0.020076_wp, -0.023889_wp, &
                      -0.023844_wp, -0.020082_wp, 0.062355_wp, 0.058806_wp, 0.041412_wp, &
                      0.049485_wp, 0.049429_wp, 0.041472_wp, 0.062267_wp, 0.058738_wp]
         case ("resp")
            charges = [ &
                      -0.60293_wp, 0.913683_wp, -0.602674_wp, -0.198556_wp, -0.063986_wp, &
                      -0.054353_wp, -0.189515_wp, 0.128751_wp, 0.131946_wp, 0.070443_wp, &
                      0.074601_wp, 0.070512_wp, 0.068296_wp, 0.12493_wp, 0.128852_wp]
         case ("mbis")
            charges = [ &
                      -0.684722_wp, 1.223113_wp, -0.684743_wp, -0.463485_wp, -0.184984_wp, &
                      -0.184615_wp, -0.463169_wp, 0.217245_wp, 0.214303_wp, 0.141753_wp, &
                      0.148298_wp, 0.148123_wp, 0.141641_wp, 0.217068_wp, 0.214161_wp]
         case ("chelpg")
            charges = [ &
                      -0.604476_wp, 0.920959_wp, -0.604202_wp, -0.206113_wp, -0.065784_wp, &
                      -0.055369_wp, -0.197225_wp, 0.130934_wp, 0.134281_wp, 0.071915_wp, &
                      0.075726_wp, 0.071446_wp, 0.069594_wp, 0.12712_wp, 0.131193_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.388775_wp, 0.438639_wp, -0.388749_wp, -0.019822_wp, -0.023685_wp, &
                      -0.023638_wp, -0.019828_wp, 0.062519_wp, 0.059067_wp, 0.041574_wp, &
                      0.049838_wp, 0.04978_wp, 0.04164_wp, 0.062433_wp, 0.058994_wp]
         case ("resp")
            charges = [ &
                      -0.605057_wp, 0.913544_wp, -0.604801_wp, -0.197068_wp, -0.063423_wp, &
                      -0.053802_wp, -0.18804_wp, 0.128611_wp, 0.132048_wp, 0.07041_wp, &
                      0.074823_wp, 0.070734_wp, 0.068272_wp, 0.1248_wp, 0.128949_wp]
         case ("mbis")
            charges = [ &
                      -0.686576_wp, 1.223897_wp, -0.686598_wp, -0.463447_wp, -0.184975_wp, &
                      -0.184609_wp, -0.463127_wp, 0.217495_wp, 0.214722_wp, 0.141955_wp, &
                      0.148842_wp, 0.148663_wp, 0.141853_wp, 0.217319_wp, 0.214573_wp]
         case ("chelpg")
            charges = [ &
                      -0.606602_wp, 0.920817_wp, -0.606329_wp, -0.204624_wp, -0.065212_wp, &
                      -0.054805_wp, -0.195749_wp, 0.130793_wp, 0.134382_wp, 0.071878_wp, &
                      0.075944_wp, 0.071664_wp, 0.069566_wp, 0.126989_wp, 0.131289_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (166)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.001088_wp, -0.001082_wp, -0.037589_wp, -0.037595_wp, -0.043537_wp, &
                      -0.043539_wp, -0.036146_wp, -0.043211_wp, -0.043214_wp, -0.036138_wp, &
                      0.026816_wp, 0.024836_wp, 0.026812_wp, 0.024835_wp, 0.030858_wp, &
                      0.030855_wp, 0.021609_wp, 0.024144_wp, 0.033315_wp, 0.033314_wp, &
                      0.024143_wp, 0.021604_wp]
         case ("resp")
            charges = [ &
                      0.099892_wp, 0.10412_wp, -0.126771_wp, -0.123808_wp, -0.242458_wp, &
                      -0.247819_wp, -0.022452_wp, -0.15401_wp, -0.136107_wp, -0.021054_wp, &
                      0.052053_wp, 0.057936_wp, 0.052079_wp, 0.05721_wp, 0.140881_wp, &
                      0.141827_wp, 0.032399_wp, 0.015009_wp, 0.139399_wp, 0.134476_wp, &
                      0.015113_wp, 0.032083_wp]
         case ("mbis")
            charges = [ &
                      0.052486_wp, 0.052504_wp, -0.239235_wp, -0.239213_wp, -0.190583_wp, &
                      -0.190549_wp, -0.203741_wp, -0.140789_wp, -0.140799_wp, -0.203642_wp, &
                      0.118333_wp, 0.116334_wp, 0.118323_wp, 0.116334_wp, 0.135224_wp, &
                      0.135208_wp, 0.110739_wp, 0.107284_wp, 0.133923_wp, 0.133922_wp, &
                      0.107268_wp, 0.110672_wp]
         case ("chelpg")
            charges = [ &
                      0.107555_wp, 0.111774_wp, -0.140003_wp, -0.137033_wp, -0.248602_wp, &
                      -0.254214_wp, -0.021938_wp, -0.153124_wp, -0.134665_wp, -0.020455_wp, &
                      0.05555_wp, 0.061187_wp, 0.055598_wp, 0.060463_wp, 0.142662_wp, &
                      0.143635_wp, 0.033187_wp, 0.015749_wp, 0.139523_wp, 0.134466_wp, &
                      0.015828_wp, 0.032856_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.004106_wp, -0.004098_wp, -0.037928_wp, -0.037936_wp, -0.048615_wp, &
                      -0.048618_wp, -0.036343_wp, -0.049754_wp, -0.049757_wp, -0.036335_wp, &
                      0.029321_wp, 0.027072_wp, 0.029318_wp, 0.027071_wp, 0.035063_wp, &
                      0.03506_wp, 0.022391_wp, 0.025977_wp, 0.036927_wp, 0.036925_wp, &
                      0.025976_wp, 0.022393_wp]
         case ("resp")
            charges = [ &
                      0.09856_wp, 0.102834_wp, -0.134518_wp, -0.131583_wp, -0.259222_wp, &
                      -0.264504_wp, -0.022298_wp, -0.169568_wp, -0.151559_wp, -0.020663_wp, &
                      0.058106_wp, 0.064247_wp, 0.058146_wp, 0.063549_wp, 0.153917_wp, &
                      0.154822_wp, 0.033144_wp, 0.017924_wp, 0.151439_wp, 0.146486_wp, &
                      0.017965_wp, 0.032776_wp]
         case ("mbis")
            charges = [ &
                      0.050737_wp, 0.050755_wp, -0.244259_wp, -0.244246_wp, -0.204844_wp, &
                      -0.204817_wp, -0.206472_wp, -0.155192_wp, -0.155201_wp, -0.206381_wp, &
                      0.123737_wp, 0.121213_wp, 0.123731_wp, 0.121214_wp, 0.146772_wp, &
                      0.146757_wp, 0.112492_wp, 0.110934_wp, 0.144859_wp, 0.144857_wp, &
                      0.110917_wp, 0.112439_wp]
         case ("chelpg")
            charges = [ &
                      0.106317_wp, 0.110674_wp, -0.148101_wp, -0.145143_wp, -0.265479_wp, &
                      -0.270907_wp, -0.02177_wp, -0.168641_wp, -0.150196_wp, -0.020031_wp, &
                      0.061699_wp, 0.067588_wp, 0.061762_wp, 0.066892_wp, 0.155725_wp, &
                      0.156638_wp, 0.033956_wp, 0.01868_wp, 0.151568_wp, 0.146507_wp, &
                      0.018692_wp, 0.033569_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.006164_wp, -0.006155_wp, -0.038089_wp, -0.038098_wp, -0.051959_wp, &
                      -0.051964_wp, -0.036353_wp, -0.05424_wp, -0.054244_wp, -0.036344_wp, &
                      0.031037_wp, 0.028599_wp, 0.031034_wp, 0.028598_wp, 0.03785_wp, &
                      0.037847_wp, 0.022808_wp, 0.0274_wp, 0.039113_wp, 0.039111_wp, &
                      0.027398_wp, 0.022816_wp]
         case ("resp")
            charges = [ &
                      0.096635_wp, 0.100913_wp, -0.139153_wp, -0.136251_wp, -0.269476_wp, &
                      -0.274722_wp, -0.022039_wp, -0.180127_wp, -0.162036_wp, -0.020231_wp, &
                      0.062213_wp, 0.068494_wp, 0.062264_wp, 0.067822_wp, 0.162385_wp, &
                      0.163265_wp, 0.033312_wp, 0.020299_wp, 0.159104_wp, 0.154131_wp, &
                      0.020292_wp, 0.032908_wp]
         case ("mbis")
            charges = [ &
                      0.049305_wp, 0.049324_wp, -0.247472_wp, -0.247465_wp, -0.213968_wp, &
                      -0.213945_wp, -0.208015_wp, -0.164909_wp, -0.164918_wp, -0.207932_wp, &
                      0.127384_wp, 0.12449_wp, 0.127379_wp, 0.124492_wp, 0.15435_wp, &
                      0.154336_wp, 0.113372_wp, 0.11365_wp, 0.151792_wp, 0.15179_wp, &
                      0.113632_wp, 0.113331_wp]
         case ("chelpg")
            charges = [ &
                      0.104445_wp, 0.108863_wp, -0.152943_wp, -0.150003_wp, -0.275785_wp, &
                      -0.281116_wp, -0.0215_wp, -0.179174_wp, -0.160715_wp, -0.01957_wp, &
                      0.065859_wp, 0.071886_wp, 0.065933_wp, 0.071214_wp, 0.164201_wp, &
                      0.165078_wp, 0.034137_wp, 0.021062_wp, 0.159232_wp, 0.154164_wp, &
                      0.021022_wp, 0.03371_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (167)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.035073_wp, -0.053467_wp, -0.055012_wp, -0.055009_wp, -0.053466_wp, &
                      0.04935_wp, 0.041591_wp, 0.041592_wp, 0.049351_wp]
         case ("resp")
            charges = [ &
                      0.05452_wp, -0.234568_wp, -0.134163_wp, -0.134854_wp, -0.234178_wp, &
                      0.198249_wp, 0.143362_wp, 0.143436_wp, 0.198196_wp]
         case ("mbis")
            charges = [ &
                      0.100247_wp, -0.256599_wp, -0.13621_wp, -0.136206_wp, -0.256615_wp, &
                      0.192678_wp, 0.15001_wp, 0.150014_wp, 0.192683_wp]
         case ("chelpg")
            charges = [ &
                      0.05521_wp, -0.235754_wp, -0.133773_wp, -0.134479_wp, -0.235355_wp, &
                      0.198703_wp, 0.143362_wp, 0.14344_wp, 0.198647_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.026828_wp, -0.056347_wp, -0.059221_wp, -0.059219_wp, -0.056347_wp, &
                      0.05574_wp, 0.046412_wp, 0.046414_wp, 0.055741_wp]
         case ("resp")
            charges = [ &
                      0.051957_wp, -0.251488_wp, -0.145145_wp, -0.14583_wp, -0.251116_wp, &
                      0.214915_wp, 0.15588_wp, 0.155954_wp, 0.214873_wp]
         case ("mbis")
            charges = [ &
                      0.092521_wp, -0.266589_wp, -0.148049_wp, -0.148044_wp, -0.266608_wp, &
                      0.206489_wp, 0.161891_wp, 0.161896_wp, 0.206494_wp]
         case ("chelpg")
            charges = [ &
                      0.052638_wp, -0.252657_wp, -0.144773_wp, -0.145472_wp, -0.252276_wp, &
                      0.215365_wp, 0.155889_wp, 0.155966_wp, 0.21532_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.02066_wp, -0.058217_wp, -0.06172_wp, -0.061717_wp, -0.058217_wp, &
                      0.059902_wp, 0.049703_wp, 0.049704_wp, 0.059903_wp]
         case ("resp")
            charges = [ &
                      0.048801_wp, -0.26219_wp, -0.151966_wp, -0.152653_wp, -0.261828_wp, &
                      0.22571_wp, 0.164189_wp, 0.164263_wp, 0.225675_wp]
         case ("mbis")
            charges = [ &
                      0.086179_wp, -0.272861_wp, -0.155377_wp, -0.155371_wp, -0.272882_wp, &
                      0.215416_wp, 0.169736_wp, 0.169741_wp, 0.215421_wp]
         case ("chelpg")
            charges = [ &
                      0.049478_wp, -0.263343_wp, -0.151611_wp, -0.15231_wp, -0.262973_wp, &
                      0.226156_wp, 0.164204_wp, 0.164282_wp, 0.226117_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (168)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.063081_wp, 0.003586_wp, -0.042834_wp, -0.03007_wp, -0.041128_wp, &
                      -0.029001_wp, -0.040893_wp, 0.054849_wp, 0.035683_wp, 0.03829_wp, &
                      0.036432_wp, 0.03884_wp, 0.03933_wp]
         case ("resp")
            charges = [ &
                      -0.25936_wp, 0.035844_wp, -0.141855_wp, -0.152949_wp, -0.051047_wp, &
                      -0.229333_wp, -0.062108_wp, 0.188919_wp, 0.129599_wp, 0.130425_wp, &
                      0.10845_wp, 0.143033_wp, 0.160382_wp]
         case ("mbis")
            charges = [ &
                      -0.228273_wp, 0.061555_wp, -0.16383_wp, -0.118194_wp, -0.155826_wp, &
                      -0.127938_wp, -0.143008_wp, 0.171129_wp, 0.140071_wp, 0.135644_wp, &
                      0.138537_wp, 0.137598_wp, 0.152538_wp]
         case ("chelpg")
            charges = [ &
                      -0.260661_wp, 0.04211_wp, -0.147659_wp, -0.150113_wp, -0.05253_wp, &
                      -0.227738_wp, -0.066572_wp, 0.189061_wp, 0.13088_wp, 0.130174_wp, &
                      0.10856_wp, 0.143116_wp, 0.161372_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.081192_wp, -0.00357_wp, -0.045526_wp, -0.03433_wp, -0.04542_wp, &
                      -0.033755_wp, -0.046008_wp, 0.069338_wp, 0.04466_wp, 0.044428_wp, &
                      0.041929_wp, 0.044694_wp, 0.044753_wp]
         case ("resp")
            charges = [ &
                      -0.301828_wp, 0.034601_wp, -0.155008_wp, -0.168705_wp, -0.06292_wp, &
                      -0.244813_wp, -0.077648_wp, 0.223218_wp, 0.14969_wp, 0.14601_wp, &
                      0.12254_wp, 0.158578_wp, 0.176286_wp]
         case ("mbis")
            charges = [ &
                      -0.267_wp, 0.059461_wp, -0.176783_wp, -0.131272_wp, -0.168219_wp, &
                      -0.140605_wp, -0.157811_wp, 0.20249_wp, 0.159243_wp, 0.15015_wp, &
                      0.151999_wp, 0.151703_wp, 0.166649_wp]
         case ("chelpg")
            charges = [ &
                      -0.303236_wp, 0.041843_wp, -0.161796_wp, -0.164842_wp, -0.065548_wp, &
                      -0.242017_wp, -0.083294_wp, 0.223356_wp, 0.151136_wp, 0.145583_wp, &
                      0.122874_wp, 0.158434_wp, 0.177507_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.088173_wp, -0.006281_wp, -0.04617_wp, -0.035562_wp, -0.046796_wp, &
                      -0.035504_wp, -0.04823_wp, 0.074931_wp, 0.048168_wp, 0.046738_wp, &
                      0.04392_wp, 0.046702_wp, 0.04626_wp]
         case ("resp")
            charges = [ &
                      -0.318399_wp, 0.034007_wp, -0.158528_wp, -0.173992_wp, -0.067398_wp, &
                      -0.249701_wp, -0.084241_wp, 0.236653_wp, 0.156974_wp, 0.151655_wp, &
                      0.127658_wp, 0.163962_wp, 0.181351_wp]
         case ("mbis")
            charges = [ &
                      -0.281727_wp, 0.05845_wp, -0.180759_wp, -0.135545_wp, -0.172378_wp, &
                      -0.144971_wp, -0.163691_wp, 0.214407_wp, 0.166405_wp, 0.155361_wp, &
                      0.156747_wp, 0.156575_wp, 0.17113_wp]
         case ("chelpg")
            charges = [ &
                      -0.319842_wp, 0.041572_wp, -0.165638_wp, -0.169787_wp, -0.070414_wp, &
                      -0.246496_wp, -0.090286_wp, 0.23679_wp, 0.158474_wp, 0.151169_wp, &
                      0.128068_wp, 0.163739_wp, 0.182651_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (169)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.008926_wp, -0.043336_wp, -0.036899_wp, -0.042145_wp, -0.03663_wp, &
                      -0.04244_wp, 0.031187_wp, 0.034987_wp, 0.03429_wp, 0.034932_wp, &
                      0.032557_wp, -0.067694_wp, 0.031022_wp, 0.031024_wp, 0.030216_wp]
         case ("resp")
            charges = [ &
                      0.41037_wp, -0.320855_wp, -0.078461_wp, -0.194966_wp, -0.083037_wp, &
                      -0.315487_wp, 0.154938_wp, 0.125479_wp, 0.137227_wp, 0.126255_wp, &
                      0.156793_wp, -0.545447_wp, 0.142398_wp, 0.142271_wp, 0.142522_wp]
         case ("mbis")
            charges = [ &
                      0.153179_wp, -0.201673_wp, -0.121159_wp, -0.156362_wp, -0.12159_wp, &
                      -0.202743_wp, 0.135074_wp, 0.132803_wp, 0.136221_wp, 0.132788_wp, &
                      0.140001_wp, -0.463955_wp, 0.146357_wp, 0.146351_wp, 0.144704_wp]
         case ("chelpg")
            charges = [ &
                      0.439472_wp, -0.340012_wp, -0.062983_wp, -0.211222_wp, -0.065684_wp, &
                      -0.338656_wp, 0.158211_wp, 0.123129_wp, 0.140207_wp, 0.123517_wp, &
                      0.161757_wp, -0.57227_wp, 0.148088_wp, 0.147996_wp, 0.148449_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.006173_wp, -0.04759_wp, -0.041685_wp, -0.04733_wp, -0.041314_wp, &
                      -0.046177_wp, 0.035382_wp, 0.038754_wp, 0.037819_wp, 0.03874_wp, &
                      0.036723_wp, -0.06849_wp, 0.033142_wp, 0.033148_wp, 0.032703_wp]
         case ("resp")
            charges = [ &
                      0.413243_wp, -0.334667_wp, -0.091483_wp, -0.206857_wp, -0.096162_wp, &
                      -0.327801_wp, 0.167063_wp, 0.136587_wp, 0.147775_wp, 0.137373_wp, &
                      0.168692_wp, -0.563764_wp, 0.149981_wp, 0.149858_wp, 0.150163_wp]
         case ("mbis")
            charges = [ &
                      0.155194_wp, -0.215405_wp, -0.131764_wp, -0.168766_wp, -0.1325_wp, &
                      -0.214546_wp, 0.146178_wp, 0.142898_wp, 0.146201_wp, 0.142937_wp, &
                      0.150531_wp, -0.475069_wp, 0.15167_wp, 0.151672_wp, 0.150767_wp]
         case ("chelpg")
            charges = [ &
                      0.441921_wp, -0.353219_wp, -0.076726_wp, -0.222406_wp, -0.079498_wp, &
                      -0.350406_wp, 0.170226_wp, 0.134388_wp, 0.15062_wp, 0.134776_wp, &
                      0.17356_wp, -0.590706_wp, 0.155714_wp, 0.155626_wp, 0.156131_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.003849_wp, -0.051124_wp, -0.045812_wp, -0.051871_wp, -0.045308_wp, &
                      -0.049229_wp, 0.038971_wp, 0.041819_wp, 0.040624_wp, 0.041868_wp, &
                      0.040323_wp, -0.069_wp, 0.034996_wp, 0.035006_wp, 0.034885_wp]
         case ("resp")
            charges = [ &
                      0.412809_wp, -0.344831_wp, -0.103057_wp, -0.216737_wp, -0.10774_wp, &
                      -0.336626_wp, 0.17709_wp, 0.145859_wp, 0.156407_wp, 0.14667_wp, &
                      0.178569_wp, -0.576908_wp, 0.156137_wp, 0.156018_wp, 0.156337_wp]
         case ("mbis")
            charges = [ &
                      0.156367_wp, -0.226654_wp, -0.140764_wp, -0.17946_wp, -0.141673_wp, &
                      -0.224104_wp, 0.155534_wp, 0.151224_wp, 0.154376_wp, 0.151327_wp, &
                      0.159446_wp, -0.483703_wp, 0.156102_wp, 0.156111_wp, 0.155868_wp]
         case ("chelpg")
            charges = [ &
                      0.441145_wp, -0.362916_wp, -0.088851_wp, -0.231746_wp, -0.091604_wp, &
                      -0.358794_wp, 0.180171_wp, 0.143773_wp, 0.159151_wp, 0.14418_wp, &
                      0.183365_wp, -0.603922_wp, 0.161899_wp, 0.161815_wp, 0.162334_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (170)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.015898_wp, -0.011111_wp, -0.039982_wp, -0.038848_wp, -0.038844_wp, &
                      -0.039995_wp, -0.01111_wp, 0.015898_wp, -0.039981_wp, -0.038848_wp, &
                      -0.038843_wp, -0.039992_wp, 0.019956_wp, 0.017556_wp, 0.020978_wp, &
                      0.01798_wp, 0.017966_wp, 0.020944_wp, 0.019952_wp, 0.01753_wp, &
                      0.019957_wp, 0.017557_wp, 0.017973_wp, 0.020971_wp, 0.017968_wp, &
                      0.020947_wp, 0.019953_wp, 0.017531_wp]
         case ("resp")
            charges = [ &
                      -0.019493_wp, 0.19722_wp, -0.159431_wp, -0.011653_wp, 0.009376_wp, &
                      -0.163108_wp, 0.198012_wp, -0.020798_wp, -0.160089_wp, -0.011889_wp, &
                      -4.6e-05_wp, -0.170824_wp, 0.029834_wp, 0.02956_wp, 0.011557_wp, &
                      0.009319_wp, 0.003553_wp, 0.00665_wp, 0.031082_wp, 0.027253_wp, &
                      0.030368_wp, 0.029165_wp, 0.01058_wp, 0.012087_wp, 0.007479_wp, &
                      0.01019_wp, 0.033655_wp, 0.03039_wp]
         case ("mbis")
            charges = [ &
                      0.094033_wp, -0.043141_wp, -0.240189_wp, -0.203127_wp, -0.202867_wp, &
                      -0.240379_wp, -0.043148_wp, 0.094035_wp, -0.240247_wp, -0.203045_wp, &
                      -0.202864_wp, -0.240398_wp, 0.106992_wp, 0.105401_wp, 0.103972_wp, &
                      0.101576_wp, 0.101475_wp, 0.103856_wp, 0.106989_wp, 0.105405_wp, &
                      0.107005_wp, 0.105415_wp, 0.101535_wp, 0.103934_wp, 0.101479_wp, &
                      0.103859_wp, 0.106994_wp, 0.105411_wp]
         case ("chelpg")
            charges = [ &
                      -0.021758_wp, 0.208408_wp, -0.169889_wp, -0.011838_wp, 0.013404_wp, &
                      -0.175519_wp, 0.209629_wp, -0.023175_wp, -0.170867_wp, -0.011727_wp, &
                      0.002774_wp, -0.183333_wp, 0.032523_wp, 0.031685_wp, 0.012019_wp, &
                      0.009695_wp, 0.003074_wp, 0.006288_wp, 0.033923_wp, 0.029603_wp, &
                      0.033086_wp, 0.031395_wp, 0.01098_wp, 0.012544_wp, 0.007364_wp, &
                      0.010177_wp, 0.03664_wp, 0.032895_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.016734_wp, -0.010924_wp, -0.040596_wp, -0.039526_wp, -0.039523_wp, &
                      -0.040609_wp, -0.010923_wp, 0.016734_wp, -0.040595_wp, -0.039525_wp, &
                      -0.039523_wp, -0.040606_wp, 0.020028_wp, 0.017869_wp, 0.020994_wp, &
                      0.018359_wp, 0.018342_wp, 0.020962_wp, 0.020025_wp, 0.017844_wp, &
                      0.02003_wp, 0.017871_wp, 0.018353_wp, 0.020986_wp, 0.018345_wp, &
                      0.020964_wp, 0.020027_wp, 0.017845_wp]
         case ("resp")
            charges = [ &
                      -0.018827_wp, 0.207104_wp, -0.168471_wp, -0.010635_wp, 0.010515_wp, &
                      -0.172093_wp, 0.207707_wp, -0.02012_wp, -0.168895_wp, -0.011023_wp, &
                      0.001065_wp, -0.179841_wp, 0.031028_wp, 0.031313_wp, 0.010942_wp, &
                      0.009735_wp, 0.003933_wp, 0.005995_wp, 0.032272_wp, 0.028982_wp, &
                      0.031521_wp, 0.03087_wp, 0.011028_wp, 0.011495_wp, 0.007871_wp, &
                      0.009567_wp, 0.034838_wp, 0.032123_wp]
         case ("mbis")
            charges = [ &
                      0.095521_wp, -0.041631_wp, -0.243552_wp, -0.205149_wp, -0.204889_wp, &
                      -0.243747_wp, -0.041637_wp, 0.095522_wp, -0.243612_wp, -0.205062_wp, &
                      -0.204891_wp, -0.243764_wp, 0.1078_wp, 0.106649_wp, 0.104553_wp, &
                      0.102826_wp, 0.102719_wp, 0.104439_wp, 0.107802_wp, 0.106655_wp, &
                      0.107814_wp, 0.106663_wp, 0.102785_wp, 0.104513_wp, 0.102723_wp, &
                      0.104443_wp, 0.107808_wp, 0.10666_wp]
         case ("chelpg")
            charges = [ &
                      -0.021088_wp, 0.2184_wp, -0.179154_wp, -0.010739_wp, 0.014651_wp, &
                      -0.184726_wp, 0.21943_wp, -0.022494_wp, -0.179887_wp, -0.010791_wp, &
                      0.003985_wp, -0.192568_wp, 0.033779_wp, 0.033491_wp, 0.011386_wp, &
                      0.010096_wp, 0.003433_wp, 0.005608_wp, 0.035172_wp, 0.031383_wp, &
                      0.0343_wp, 0.033152_wp, 0.011416_wp, 0.011938_wp, 0.007737_wp, &
                      0.009532_wp, 0.037881_wp, 0.034678_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.01765_wp, -0.010686_wp, -0.041132_wp, -0.040132_wp, -0.04013_wp, &
                      -0.041144_wp, -0.010685_wp, 0.01765_wp, -0.04113_wp, -0.04013_wp, &
                      -0.040131_wp, -0.041141_wp, 0.020006_wp, 0.018179_wp, 0.020883_wp, &
                      0.018746_wp, 0.018725_wp, 0.020852_wp, 0.020006_wp, 0.018155_wp, &
                      0.020009_wp, 0.018181_wp, 0.01874_wp, 0.020875_wp, 0.018728_wp, &
                      0.020855_wp, 0.020008_wp, 0.018155_wp]
         case ("resp")
            charges = [ &
                      -0.018111_wp, 0.216906_wp, -0.176909_wp, -0.009448_wp, 0.011815_wp, &
                      -0.180482_wp, 0.217325_wp, -0.019395_wp, -0.177104_wp, -0.009981_wp, &
                      0.002337_wp, -0.188255_wp, 0.031916_wp, 0.032915_wp, 0.010057_wp, &
                      0.010125_wp, 0.004289_wp, 0.005073_wp, 0.033159_wp, 0.030562_wp, &
                      0.03237_wp, 0.032424_wp, 0.011449_wp, 0.010632_wp, 0.008237_wp, &
                      0.008675_wp, 0.035716_wp, 0.033705_wp]
         case ("mbis")
            charges = [ &
                      0.097153_wp, -0.040255_wp, -0.246518_wp, -0.206891_wp, -0.206629_wp, &
                      -0.246719_wp, -0.040259_wp, 0.097153_wp, -0.24658_wp, -0.2068_wp, &
                      -0.206635_wp, -0.246734_wp, 0.108352_wp, 0.107816_wp, 0.104832_wp, &
                      0.104033_wp, 0.103919_wp, 0.104721_wp, 0.108358_wp, 0.107824_wp, &
                      0.108367_wp, 0.107831_wp, 0.103991_wp, 0.104791_wp, 0.103924_wp, &
                      0.104726_wp, 0.108364_wp, 0.107828_wp]
         case ("chelpg")
            charges = [ &
                      -0.020368_wp, 0.22829_wp, -0.187788_wp, -0.009465_wp, 0.016063_wp, &
                      -0.19331_wp, 0.229132_wp, -0.021764_wp, -0.18828_wp, -0.009675_wp, &
                      0.005365_wp, -0.201176_wp, 0.034721_wp, 0.035138_wp, 0.010478_wp, &
                      0.010467_wp, 0.003765_wp, 0.004658_wp, 0.036111_wp, 0.033005_wp, &
                      0.0352_wp, 0.034748_wp, 0.011821_wp, 0.011056_wp, 0.008081_wp, &
                      0.008614_wp, 0.038811_wp, 0.036303_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (171)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070428_wp, -0.033069_wp, -0.034636_wp, 0.057181_wp, -0.179526_wp, &
                      0.477353_wp, -0.414916_wp, -0.182872_wp, 0.065768_wp, -0.036487_wp, &
                      -0.03169_wp, -0.0734_wp, -0.188203_wp, 0.054802_wp, -0.036762_wp, &
                      -0.030924_wp, -0.070138_wp, 0.019471_wp, 0.027577_wp, 0.023893_wp, &
                      0.023288_wp, 0.028534_wp, 0.02727_wp, 0.030416_wp, 0.028202_wp, &
                      0.029141_wp, 0.038195_wp, 0.04216_wp, 0.023384_wp, 0.032051_wp, &
                      0.023506_wp, 0.024252_wp, 0.018943_wp, 0.024799_wp, 0.020709_wp, &
                      0.024731_wp, 0.032906_wp, 0.034036_wp, 0.025798_wp, 0.023385_wp, &
                      0.027628_wp, 0.021482_wp, 0.028464_wp, 0.02377_wp]
         case ("resp")
            charges = [ &
                      -0.283634_wp, 0.205623_wp, -0.227214_wp, 0.211163_wp, -0.309923_wp, &
                      0.86639_wp, -0.609421_wp, -0.255589_wp, 0.004484_wp, -0.117419_wp, &
                      0.07152_wp, -0.263925_wp, -0.390339_wp, 0.234026_wp, -0.19106_wp, &
                      0.03128_wp, -0.222055_wp, 0.061213_wp, 0.066466_wp, 0.071116_wp, &
                      -0.019028_wp, -0.01035_wp, 0.064916_wp, 0.065541_wp, -0.002018_wp, &
                      0.001809_wp, 0.079143_wp, 0.089482_wp, 0.087997_wp, 0.038394_wp, &
                      0.004734_wp, 0.022054_wp, 0.062473_wp, 0.062242_wp, 0.083443_wp, &
                      0.000613_wp, 0.057276_wp, 0.070818_wp, 0.075504_wp, 0.011266_wp, &
                      0.035442_wp, 0.046283_wp, 0.069094_wp, 0.050168_wp]
         case ("mbis")
            charges = [ &
                      -0.422466_wp, -0.150177_wp, -0.290397_wp, 0.105058_wp, -0.567292_wp, &
                      1.759358_wp, -0.854415_wp, -0.606174_wp, 0.11699_wp, -0.322129_wp, &
                      -0.12875_wp, -0.422713_wp, -0.589201_wp, 0.119964_wp, -0.304442_wp, &
                      -0.157629_wp, -0.421996_wp, 0.121567_wp, 0.134644_wp, 0.13671_wp, &
                      0.100545_wp, 0.12515_wp, 0.120926_wp, 0.130583_wp, 0.095875_wp, &
                      0.074556_wp, 0.095588_wp, 0.100646_wp, 0.134646_wp, 0.130788_wp, &
                      0.099662_wp, 0.100457_wp, 0.11866_wp, 0.129234_wp, 0.145489_wp, &
                      0.081047_wp, 0.080396_wp, 0.14321_wp, 0.129555_wp, 0.103662_wp, &
                      0.112485_wp, 0.12378_wp, 0.13752_wp, 0.129071_wp]
         case ("chelpg")
            charges = [ &
                      -0.304307_wp, 0.225175_wp, -0.247528_wp, 0.235726_wp, -0.321717_wp, &
                      0.889416_wp, -0.613541_wp, -0.268394_wp, 0.014218_wp, -0.13746_wp, &
                      0.092427_wp, -0.285036_wp, -0.402494_wp, 0.261983_wp, -0.215675_wp, &
                      0.039369_wp, -0.235711_wp, 0.066007_wp, 0.070671_wp, 0.075652_wp, &
                      -0.022077_wp, -0.013434_wp, 0.068031_wp, 0.069191_wp, -0.007438_wp, &
                      -0.004334_wp, 0.077048_wp, 0.088836_wp, 0.093202_wp, 0.041726_wp, &
                      0.001392_wp, 0.018545_wp, 0.067587_wp, 0.065903_wp, 0.090841_wp, &
                      -0.005282_wp, 0.052035_wp, 0.075875_wp, 0.082097_wp, 0.010729_wp, &
                      0.035588_wp, 0.049047_wp, 0.072466_wp, 0.053647_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.072368_wp, -0.035006_wp, -0.034159_wp, 0.061944_wp, -0.191502_wp, &
                      0.469046_wp, -0.46639_wp, -0.180354_wp, 0.066352_wp, -0.033291_wp, &
                      -0.032004_wp, -0.075629_wp, -0.200419_wp, 0.059595_wp, -0.037023_wp, &
                      -0.029426_wp, -0.070826_wp, 0.023799_wp, 0.026721_wp, 0.022872_wp, &
                      0.025902_wp, 0.020383_wp, 0.03596_wp, 0.030895_wp, 0.034819_wp, &
                      0.041518_wp, 0.040827_wp, 0.043905_wp, 0.02989_wp, 0.038453_wp, &
                      0.025496_wp, 0.024808_wp, 0.022427_wp, 0.023845_wp, 0.018929_wp, &
                      0.034197_wp, 0.042417_wp, 0.030872_wp, 0.029749_wp, 0.025787_wp, &
                      0.028415_wp, 0.025452_wp, 0.027309_wp, 0.025853_wp]
         case ("resp")
            charges = [ &
                      -0.295891_wp, 0.214988_wp, -0.237769_wp, 0.235436_wp, -0.370158_wp, &
                      0.957577_wp, -0.707222_wp, -0.279345_wp, 0.018362_wp, -0.097613_wp, &
                      0.061504_wp, -0.278749_wp, -0.44779_wp, 0.252867_wp, -0.197931_wp, &
                      0.051155_wp, -0.236896_wp, 0.071232_wp, 0.066378_wp, 0.067831_wp, &
                      -0.014866_wp, -0.027016_wp, 0.081243_wp, 0.065754_wp, 0.010309_wp, &
                      0.019801_wp, 0.082175_wp, 0.087533_wp, 0.095686_wp, 0.042909_wp, &
                      0.008648_wp, 0.024372_wp, 0.071956_wp, 0.063073_wp, 0.083718_wp, &
                      0.020153_wp, 0.072211_wp, 0.063332_wp, 0.083556_wp, 0.010023_wp, &
                      0.031949_wp, 0.054577_wp, 0.067004_wp, 0.053936_wp]
         case ("mbis")
            charges = [ &
                      -0.435503_wp, -0.14678_wp, -0.301418_wp, 0.106867_wp, -0.593192_wp, &
                      1.782174_wp, -0.924526_wp, -0.605888_wp, 0.112037_wp, -0.320477_wp, &
                      -0.134033_wp, -0.433144_wp, -0.614896_wp, 0.120197_wp, -0.311494_wp, &
                      -0.152152_wp, -0.434524_wp, 0.131911_wp, 0.135542_wp, 0.136877_wp, &
                      0.106123_wp, 0.110724_wp, 0.137571_wp, 0.132115_wp, 0.107416_wp, &
                      0.096167_wp, 0.10126_wp, 0.104231_wp, 0.144961_wp, 0.140538_wp, &
                      0.103852_wp, 0.102148_wp, 0.127162_wp, 0.129477_wp, 0.143706_wp, &
                      0.098433_wp, 0.096782_wp, 0.138281_wp, 0.138005_wp, 0.106866_wp, &
                      0.11292_wp, 0.132487_wp, 0.136963_wp, 0.134274_wp]
         case ("chelpg")
            charges = [ &
                      -0.31687_wp, 0.235025_wp, -0.258497_wp, 0.259313_wp, -0.381629_wp, &
                      0.980244_wp, -0.71131_wp, -0.292394_wp, 0.029561_wp, -0.117663_wp, &
                      0.08175_wp, -0.299751_wp, -0.459541_wp, 0.279935_wp, -0.222969_wp, &
                      0.06072_wp, -0.251484_wp, 0.07608_wp, 0.070627_wp, 0.072491_wp, &
                      -0.018006_wp, -0.0302_wp, 0.084497_wp, 0.06953_wp, 0.005132_wp, &
                      0.013909_wp, 0.079703_wp, 0.08648_wp, 0.100899_wp, 0.046129_wp, &
                      0.005521_wp, 0.02102_wp, 0.077065_wp, 0.066757_wp, 0.091052_wp, &
                      0.014415_wp, 0.067194_wp, 0.068484_wp, 0.090203_wp, 0.009165_wp, &
                      0.03173_wp, 0.057557_wp, 0.070532_wp, 0.057597_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.07265_wp, -0.035413_wp, -0.033998_wp, 0.062885_wp, -0.193851_wp, &
                      0.467056_wp, -0.475276_wp, -0.179766_wp, 0.066399_wp, -0.032605_wp, &
                      -0.032029_wp, -0.075921_wp, -0.202636_wp, 0.060483_wp, -0.037079_wp, &
                      -0.029132_wp, -0.070824_wp, 0.024916_wp, 0.026385_wp, 0.022521_wp, &
                      0.02639_wp, 0.018417_wp, 0.037851_wp, 0.03079_wp, 0.036233_wp, &
                      0.043934_wp, 0.041094_wp, 0.043863_wp, 0.031419_wp, 0.03958_wp, &
                      0.025919_wp, 0.024634_wp, 0.023404_wp, 0.02352_wp, 0.018487_wp, &
                      0.036234_wp, 0.044118_wp, 0.029962_wp, 0.030455_wp, 0.026209_wp, &
                      0.028381_wp, 0.026434_wp, 0.02694_wp, 0.026305_wp]
         case ("resp")
            charges = [ &
                      -0.297275_wp, 0.215737_wp, -0.238964_wp, 0.239777_wp, -0.38163_wp, &
                      0.97286_wp, -0.724206_wp, -0.283526_wp, 0.02067_wp, -0.093029_wp, &
                      0.059241_wp, -0.281266_wp, -0.458008_wp, 0.25596_wp, -0.198891_wp, &
                      0.054838_wp, -0.239008_wp, 0.073642_wp, 0.065906_wp, 0.06658_wp, &
                      -0.013782_wp, -0.030838_wp, 0.084668_wp, 0.065324_wp, 0.013285_wp, &
                      0.023415_wp, 0.082401_wp, 0.086583_wp, 0.097678_wp, 0.043426_wp, &
                      0.009626_wp, 0.024304_wp, 0.074357_wp, 0.062909_wp, 0.083578_wp, &
                      0.02464_wp, 0.074839_wp, 0.061296_wp, 0.08494_wp, 0.009782_wp, &
                      0.030939_wp, 0.056483_wp, 0.066154_wp, 0.054584_wp]
         case ("mbis")
            charges = [ &
                      -0.437526_wp, -0.145976_wp, -0.303214_wp, 0.10719_wp, -0.598017_wp, &
                      1.784593_wp, -0.936383_wp, -0.605342_wp, 0.111434_wp, -0.320141_wp, &
                      -0.134797_wp, -0.43466_wp, -0.619315_wp, 0.120015_wp, -0.312448_wp, &
                      -0.151008_wp, -0.436438_wp, 0.134343_wp, 0.135312_wp, 0.136527_wp, &
                      0.107178_wp, 0.107188_wp, 0.14109_wp, 0.131977_wp, 0.109888_wp, &
                      0.100338_wp, 0.101911_wp, 0.104235_wp, 0.147438_wp, 0.142192_wp, &
                      0.104698_wp, 0.101912_wp, 0.129271_wp, 0.129143_wp, 0.14306_wp, &
                      0.102202_wp, 0.099699_wp, 0.136739_wp, 0.139499_wp, 0.107406_wp, &
                      0.112654_wp, 0.134462_wp, 0.136458_wp, 0.135255_wp]
         case ("chelpg")
            charges = [ &
                      -0.318291_wp, 0.23585_wp, -0.259746_wp, 0.263495_wp, -0.393031_wp, &
                      0.99545_wp, -0.728286_wp, -0.296599_wp, 0.032074_wp, -0.112972_wp, &
                      0.079295_wp, -0.302207_wp, -0.469677_wp, 0.28283_wp, -0.223956_wp, &
                      0.064644_wp, -0.25374_wp, 0.078495_wp, 0.07016_wp, 0.071258_wp, &
                      -0.016939_wp, -0.034039_wp, 0.087944_wp, 0.069119_wp, 0.008158_wp, &
                      0.017579_wp, 0.079872_wp, 0.085466_wp, 0.102865_wp, 0.046605_wp, &
                      0.00655_wp, 0.020992_wp, 0.079456_wp, 0.066593_wp, 0.090889_wp, &
                      0.018939_wp, 0.069869_wp, 0.066455_wp, 0.091583_wp, 0.008869_wp, &
                      0.030658_wp, 0.059498_wp, 0.069705_wp, 0.058271_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (172)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.004831_wp, 0.059571_wp, 0.06267_wp, -0.048795_wp, -0.03935_wp, &
                      -0.029265_wp]
         case ("resp")
            charges = [ &
                      -0.111032_wp, -0.005738_wp, 0.176768_wp, -0.040207_wp, -0.0262_wp, &
                      0.006409_wp]
         case ("mbis")
            charges = [ &
                      -0.132772_wp, 0.124236_wp, 0.203881_wp, -0.070143_wp, -0.078206_wp, &
                      -0.046995_wp]
         case ("chelpg")
            charges = [ &
                      -0.113244_wp, -0.0041_wp, 0.177523_wp, -0.039753_wp, -0.026557_wp, &
                      0.006131_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.000186_wp, 0.058038_wp, 0.073752_wp, -0.055221_wp, -0.040954_wp, &
                      -0.035801_wp]
         case ("resp")
            charges = [ &
                      -0.099819_wp, -0.015586_wp, 0.19624_wp, -0.052325_wp, -0.026748_wp, &
                      -0.001763_wp]
         case ("mbis")
            charges = [ &
                      -0.128634_wp, 0.121404_wp, 0.221624_wp, -0.079679_wp, -0.079974_wp, &
                      -0.05474_wp]
         case ("chelpg")
            charges = [ &
                      -0.10173_wp, -0.014227_wp, 0.19691_wp, -0.051933_wp, -0.027043_wp, &
                      -0.001977_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.003159_wp, 0.05717_wp, 0.079399_wp, -0.05838_wp, -0.041553_wp, &
                      -0.039795_wp]
         case ("resp")
            charges = [ &
                      -0.092168_wp, -0.02106_wp, 0.205935_wp, -0.058794_wp, -0.026706_wp, &
                      -0.007206_wp]
         case ("mbis")
            charges = [ &
                      -0.125248_wp, 0.119514_wp, 0.230323_wp, -0.084581_wp, -0.080533_wp, &
                      -0.059475_wp]
         case ("chelpg")
            charges = [ &
                      -0.093924_wp, -0.019827_wp, 0.206552_wp, -0.058431_wp, -0.026975_wp, &
                      -0.007395_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (173)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.13638_wp, 0.010282_wp, 0.01033_wp, 0.01023_wp, 0.009368_wp, &
                      -0.072011_wp, -0.07197_wp, 0.009355_wp, -0.071985_wp, 0.025501_wp, &
                      0.025444_wp, 0.025429_wp, 0.009324_wp, 0.024956_wp, 0.024769_wp, &
                      0.022644_wp, 0.024965_wp, 0.02485_wp, 0.022577_wp, 0.024885_wp, &
                      0.024786_wp, 0.022647_wp]
         case ("resp")
            charges = [ &
                      -0.462184_wp, 0.239788_wp, 0.24896_wp, 0.239872_wp, -0.028745_wp, &
                      -0.433605_wp, -0.423359_wp, -0.028964_wp, -0.419177_wp, 0.042592_wp, &
                      0.041734_wp, 0.039225_wp, -0.028_wp, 0.097558_wp, 0.11207_wp, &
                      0.11999_wp, 0.09509_wp, 0.110479_wp, 0.117149_wp, 0.094032_wp, &
                      0.108677_wp, 0.116819_wp]
         case ("mbis")
            charges = [ &
                      -0.35358_wp, 0.030498_wp, 0.029912_wp, 0.030176_wp, 0.046754_wp, &
                      -0.473842_wp, -0.47377_wp, 0.046706_wp, -0.473891_wp, 0.099363_wp, &
                      0.09945_wp, 0.099375_wp, 0.046668_wp, 0.137463_wp, 0.145663_wp, &
                      0.132385_wp, 0.137426_wp, 0.145794_wp, 0.132079_wp, 0.137348_wp, &
                      0.145752_wp, 0.132264_wp]
         case ("chelpg")
            charges = [ &
                      -0.486232_wp, 0.270935_wp, 0.280921_wp, 0.271204_wp, -0.035901_wp, &
                      -0.457979_wp, -0.447204_wp, -0.036241_wp, -0.442825_wp, 0.035017_wp, &
                      0.034024_wp, 0.031421_wp, -0.035161_wp, 0.102016_wp, 0.117441_wp, &
                      0.125809_wp, 0.099477_wp, 0.115762_wp, 0.122797_wp, 0.098356_wp, &
                      0.113901_wp, 0.122461_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.143171_wp, 0.009283_wp, 0.009313_wp, 0.009204_wp, 0.011248_wp, &
                      -0.074115_wp, -0.074071_wp, 0.011232_wp, -0.074082_wp, 0.027655_wp, &
                      0.027562_wp, 0.027561_wp, 0.01121_wp, 0.026832_wp, 0.022757_wp, &
                      0.02414_wp, 0.026828_wp, 0.022827_wp, 0.024115_wp, 0.026747_wp, &
                      0.022779_wp, 0.024142_wp]
         case ("resp")
            charges = [ &
                      -0.508111_wp, 0.246057_wp, 0.255308_wp, 0.245489_wp, -0.02337_wp, &
                      -0.448087_wp, -0.437728_wp, -0.023648_wp, -0.433746_wp, 0.049001_wp, &
                      0.048223_wp, 0.045554_wp, -0.022493_wp, 0.103733_wp, 0.111841_wp, &
                      0.125892_wp, 0.101235_wp, 0.110273_wp, 0.123136_wp, 0.10018_wp, &
                      0.108511_wp, 0.122749_wp]
         case ("mbis")
            charges = [ &
                      -0.37335_wp, 0.030849_wp, 0.03021_wp, 0.030426_wp, 0.050671_wp, &
                      -0.485453_wp, -0.485371_wp, 0.05063_wp, -0.48548_wp, 0.104087_wp, &
                      0.104137_wp, 0.104074_wp, 0.050622_wp, 0.142972_wp, 0.14443_wp, &
                      0.137372_wp, 0.142906_wp, 0.144538_wp, 0.137136_wp, 0.142825_wp, &
                      0.14452_wp, 0.137241_wp]
         case ("chelpg")
            charges = [ &
                      -0.531878_wp, 0.276902_wp, 0.286956_wp, 0.27651_wp, -0.030449_wp, &
                      -0.472433_wp, -0.461555_wp, -0.030846_wp, -0.457379_wp, 0.041523_wp, &
                      0.040612_wp, 0.037849_wp, -0.029574_wp, 0.108206_wp, 0.117216_wp, &
                      0.131706_wp, 0.105638_wp, 0.115562_wp, 0.128782_wp, 0.104521_wp, &
                      0.113741_wp, 0.12839_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.149461_wp, 0.008511_wp, 0.008521_wp, 0.008403_wp, 0.013554_wp, &
                      -0.075885_wp, -0.075839_wp, 0.013538_wp, -0.075844_wp, 0.029413_wp, &
                      0.029282_wp, 0.0293_wp, 0.013527_wp, 0.028785_wp, 0.020179_wp, &
                      0.025364_wp, 0.028766_wp, 0.020234_wp, 0.025387_wp, 0.028683_wp, &
                      0.020209_wp, 0.025366_wp]
         case ("resp")
            charges = [ &
                      -0.552765_wp, 0.252238_wp, 0.261494_wp, 0.25098_wp, -0.017113_wp, &
                      -0.459515_wp, -0.449171_wp, -0.017443_wp, -0.445344_wp, 0.054124_wp, &
                      0.053443_wp, 0.05062_wp, -0.016087_wp, 0.109606_wp, 0.10986_wp, &
                      0.130893_wp, 0.107107_wp, 0.108338_wp, 0.128274_wp, 0.106041_wp, &
                      0.106616_wp, 0.127803_wp]
         case ("mbis")
            charges = [ &
                      -0.392335_wp, 0.031219_wp, 0.030518_wp, 0.030686_wp, 0.055336_wp, &
                      -0.494556_wp, -0.494469_wp, 0.055308_wp, -0.494561_wp, 0.107906_wp, &
                      0.107915_wp, 0.107876_wp, 0.055332_wp, 0.148216_wp, 0.141727_wp, &
                      0.141459_wp, 0.148123_wp, 0.141808_wp, 0.141305_wp, 0.148036_wp, &
                      0.141826_wp, 0.14132_wp]
         case ("chelpg")
            charges = [ &
                      -0.576229_wp, 0.28276_wp, 0.292806_wp, 0.281668_wp, -0.024112_wp, &
                      -0.483795_wp, -0.47294_wp, -0.024558_wp, -0.468922_wp, 0.046742_wp, &
                      0.045932_wp, 0.043014_wp, -0.023085_wp, 0.114082_wp, 0.115231_wp, &
                      0.136692_wp, 0.111516_wp, 0.113624_wp, 0.133908_wp, 0.110388_wp, &
                      0.111844_wp, 0.133433_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (174)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070361_wp, -0.034208_wp, -0.036437_wp, -0.03626_wp, -0.036207_wp, &
                      -0.036242_wp, -0.036209_wp, -0.036265_wp, -0.036442_wp, -0.034217_wp, &
                      0.023279_wp, 0.021435_wp, 0.021422_wp, 0.019059_wp, 0.019097_wp, &
                      0.018169_wp, 0.018215_wp, 0.018216_wp, 0.018194_wp, 0.018168_wp, &
                      0.018207_wp, 0.018158_wp, 0.018164_wp, 0.018204_wp, 0.018168_wp, &
                      0.018195_wp, 0.018211_wp, 0.018179_wp, 0.018196_wp, -0.070357_wp, &
                      0.019078_wp, 0.019073_wp, 0.021427_wp, 0.023279_wp, 0.02143_wp]
         case ("resp")
            charges = [ &
                      -0.329855_wp, 0.197648_wp, 0.013334_wp, -0.049094_wp, 0.065963_wp, &
                      0.028662_wp, 0.057908_wp, -0.054082_wp, 0.017621_wp, 0.206536_wp, &
                      0.07766_wp, 0.069503_wp, 0.07011_wp, -0.031818_wp, -0.033705_wp, &
                      -0.002728_wp, -0.004136_wp, -0.001536_wp, -0.001741_wp, -0.020481_wp, &
                      -0.02173_wp, -0.010864_wp, -0.010833_wp, -0.017255_wp, -0.019894_wp, &
                      0.000147_wp, -0.000604_wp, -0.003834_wp, -0.004276_wp, -0.355582_wp, &
                      -0.033913_wp, -0.033272_wp, 0.076907_wp, 0.083382_wp, 0.075851_wp]
         case ("mbis")
            charges = [ &
                      -0.42185_wp, -0.160196_wp, -0.211358_wp, -0.214007_wp, -0.208632_wp, &
                      -0.209912_wp, -0.2086_wp, -0.213982_wp, -0.211295_wp, -0.160127_wp, &
                      0.129728_wp, 0.127565_wp, 0.127532_wp, 0.099773_wp, 0.099933_wp, &
                      0.105572_wp, 0.105618_wp, 0.105314_wp, 0.105299_wp, 0.104856_wp, &
                      0.104833_wp, 0.104924_wp, 0.105061_wp, 0.104775_wp, 0.104885_wp, &
                      0.105175_wp, 0.1054_wp, 0.105619_wp, 0.105508_wp, -0.421808_wp, &
                      0.099949_wp, 0.099703_wp, 0.127469_wp, 0.129701_wp, 0.12759_wp]
         case ("chelpg")
            charges = [ &
                      -0.34713_wp, 0.210916_wp, 0.011526_wp, -0.055948_wp, 0.076596_wp, &
                      0.020939_wp, 0.068476_wp, -0.061381_wp, 0.016284_wp, 0.219696_wp, &
                      0.081635_wp, 0.07334_wp, 0.073963_wp, -0.034275_wp, -0.036166_wp, &
                      -0.002584_wp, -0.003988_wp, -0.000371_wp, -0.000596_wp, -0.022689_wp, &
                      -0.023984_wp, -0.009842_wp, -0.009811_wp, -0.019412_wp, -0.022098_wp, &
                      0.001396_wp, 0.000663_wp, -0.003769_wp, -0.004206_wp, -0.373342_wp, &
                      -0.036332_wp, -0.035698_wp, 0.080851_wp, 0.087526_wp, 0.079814_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.071713_wp, -0.034688_wp, -0.03691_wp, -0.036768_wp, -0.036744_wp, &
                      -0.036781_wp, -0.036745_wp, -0.036771_wp, -0.036914_wp, -0.034696_wp, &
                      0.023518_wp, 0.021702_wp, 0.021699_wp, 0.019421_wp, 0.019448_wp, &
                      0.018517_wp, 0.018553_wp, 0.018522_wp, 0.0185_wp, 0.01844_wp, &
                      0.018482_wp, 0.018445_wp, 0.018433_wp, 0.018484_wp, 0.018434_wp, &
                      0.018515_wp, 0.018501_wp, 0.01851_wp, 0.018549_wp, -0.071708_wp, &
                      0.01942_wp, 0.019442_wp, 0.021703_wp, 0.023522_wp, 0.021695_wp]
         case ("resp")
            charges = [ &
                      -0.340356_wp, 0.207305_wp, 0.015938_wp, -0.047879_wp, 0.071268_wp, &
                      0.032715_wp, 0.063271_wp, -0.053017_wp, 0.02066_wp, 0.216381_wp, &
                      0.07963_wp, 0.071305_wp, 0.071951_wp, -0.033905_wp, -0.035891_wp, &
                      -0.003853_wp, -0.005376_wp, -0.002947_wp, -0.003145_wp, -0.022758_wp, &
                      -0.023983_wp, -0.012798_wp, -0.012806_wp, -0.019491_wp, -0.022164_wp, &
                      -0.001221_wp, -0.002068_wp, -0.005213_wp, -0.005534_wp, -0.366543_wp, &
                      -0.036143_wp, -0.035404_wp, 0.078842_wp, 0.085469_wp, 0.077762_wp]
         case ("mbis")
            charges = [ &
                      -0.429137_wp, -0.160502_wp, -0.213357_wp, -0.216034_wp, -0.210568_wp, &
                      -0.211848_wp, -0.210533_wp, -0.216004_wp, -0.213287_wp, -0.16043_wp, &
                      0.131548_wp, 0.129426_wp, 0.129412_wp, 0.100698_wp, 0.100836_wp, &
                      0.106668_wp, 0.106691_wp, 0.106337_wp, 0.106323_wp, 0.105823_wp, &
                      0.105806_wp, 0.105919_wp, 0.106022_wp, 0.105758_wp, 0.105838_wp, &
                      0.106225_wp, 0.106393_wp, 0.106678_wp, 0.106611_wp, -0.429092_wp, &
                      0.100835_wp, 0.100641_wp, 0.129346_wp, 0.131527_wp, 0.129446_wp]
         case ("chelpg")
            charges = [ &
                      -0.357762_wp, 0.220611_wp, 0.014336_wp, -0.054834_wp, 0.082062_wp, &
                      0.025137_wp, 0.074073_wp, -0.060499_wp, 0.0196_wp, 0.22956_wp, &
                      0.083644_wp, 0.075169_wp, 0.075831_wp, -0.036376_wp, -0.038367_wp, &
                      -0.003764_wp, -0.005284_wp, -0.001788_wp, -0.002007_wp, -0.02503_wp, &
                      -0.026301_wp, -0.01184_wp, -0.011847_wp, -0.021726_wp, -0.024447_wp, &
                      3.4e-05_wp, -0.000795_wp, -0.005218_wp, -0.005535_wp, -0.38444_wp, &
                      -0.038573_wp, -0.037842_wp, 0.082812_wp, 0.089654_wp, 0.081752_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073154_wp, -0.03516_wp, -0.037392_wp, -0.037292_wp, -0.037298_wp, &
                      -0.037332_wp, -0.037297_wp, -0.037293_wp, -0.037394_wp, -0.035167_wp, &
                      0.023715_wp, 0.021956_wp, 0.021966_wp, 0.01982_wp, 0.019833_wp, &
                      0.018892_wp, 0.018913_wp, 0.018855_wp, 0.018833_wp, 0.018718_wp, &
                      0.018761_wp, 0.018744_wp, 0.018709_wp, 0.01877_wp, 0.018703_wp, &
                      0.018867_wp, 0.018815_wp, 0.018863_wp, 0.018928_wp, -0.073147_wp, &
                      0.019795_wp, 0.019851_wp, 0.021967_wp, 0.023723_wp, 0.021946_wp]
         case ("resp")
            charges = [ &
                      -0.351521_wp, 0.218411_wp, 0.019203_wp, -0.046012_wp, 0.077644_wp, &
                      0.037695_wp, 0.06971_wp, -0.051322_wp, 0.024422_wp, 0.227693_wp, &
                      0.081586_wp, 0.07311_wp, 0.073805_wp, -0.036369_wp, -0.03847_wp, &
                      -0.005295_wp, -0.006956_wp, -0.004729_wp, -0.004918_wp, -0.025536_wp, &
                      -0.026739_wp, -0.015191_wp, -0.01525_wp, -0.022224_wp, -0.024937_wp, &
                      -0.002947_wp, -0.003915_wp, -0.006948_wp, -0.00713_wp, -0.37823_wp, &
                      -0.038773_wp, -0.037915_wp, 0.080799_wp, 0.087558_wp, 0.079691_wp]
         case ("mbis")
            charges = [ &
                      -0.436818_wp, -0.160732_wp, -0.215418_wp, -0.218141_wp, -0.212562_wp, &
                      -0.213836_wp, -0.212522_wp, -0.218104_wp, -0.215342_wp, -0.160656_wp, &
                      0.133377_wp, 0.131347_wp, 0.131359_wp, 0.101684_wp, 0.101794_wp, &
                      0.107819_wp, 0.107811_wp, 0.107419_wp, 0.107406_wp, 0.106811_wp, &
                      0.106798_wp, 0.10695_wp, 0.107009_wp, 0.106763_wp, 0.10681_wp, &
                      0.107343_wp, 0.107435_wp, 0.107783_wp, 0.10777_wp, -0.436769_wp, &
                      0.101771_wp, 0.101644_wp, 0.131287_wp, 0.133363_wp, 0.131363_wp]
         case ("chelpg")
            charges = [ &
                      -0.369055_wp, 0.231734_wp, 0.017869_wp, -0.053035_wp, 0.088581_wp, &
                      0.030331_wp, 0.08074_wp, -0.058967_wp, 0.02371_wp, 0.24087_wp, &
                      0.085638_wp, 0.077_wp, 0.077711_wp, -0.038855_wp, -0.04096_wp, &
                      -0.00528_wp, -0.006937_wp, -0.003592_wp, -0.003801_wp, -0.027878_wp, &
                      -0.029127_wp, -0.014315_wp, -0.014373_wp, -0.024546_wp, -0.027306_wp, &
                      -0.001698_wp, -0.002648_wp, -0.007044_wp, -0.007222_wp, -0.396263_wp, &
                      -0.041213_wp, -0.040363_wp, 0.084797_wp, 0.091787_wp, 0.083708_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (175)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.302821_wp, 0.15141_wp, 0.15141_wp]
         case ("resp")
            charges = [ &
                      -0.702136_wp, 0.350708_wp, 0.351429_wp]
         case ("mbis")
            charges = [ &
                      -0.854498_wp, 0.427248_wp, 0.427249_wp]
         case ("chelpg")
            charges = [ &
                      -0.702403_wp, 0.350842_wp, 0.351562_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.337686_wp, 0.168843_wp, 0.168842_wp]
         case ("resp")
            charges = [ &
                      -0.825313_wp, 0.412365_wp, 0.412949_wp]
         case ("mbis")
            charges = [ &
                      -0.96398_wp, 0.481989_wp, 0.48199_wp]
         case ("chelpg")
            charges = [ &
                      -0.825583_wp, 0.4125_wp, 0.413083_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.338197_wp, 0.169099_wp, 0.169098_wp]
         case ("resp")
            charges = [ &
                      -0.827161_wp, 0.41329_wp, 0.413872_wp]
         case ("mbis")
            charges = [ &
                      -0.965562_wp, 0.48278_wp, 0.482781_wp]
         case ("chelpg")
            charges = [ &
                      -0.827427_wp, 0.413423_wp, 0.414004_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (176)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.073498_wp, -0.073507_wp, 0.028122_wp, -0.108291_wp, 0.0281_wp, &
                      0.044639_wp, 0.044648_wp, 0.054886_wp, 0.054899_wp]
         case ("resp")
            charges = [ &
                      -0.179955_wp, -0.181347_wp, -0.07318_wp, -0.121017_wp, -0.077538_wp, &
                      0.153742_wp, 0.153451_wp, 0.161923_wp, 0.163922_wp]
         case ("mbis")
            charges = [ &
                      -0.249108_wp, -0.249266_wp, 0.022429_wp, -0.186894_wp, 0.022273_wp, &
                      0.173523_wp, 0.173558_wp, 0.146721_wp, 0.146763_wp]
         case ("chelpg")
            charges = [ &
                      -0.180382_wp, -0.181897_wp, -0.072973_wp, -0.12106_wp, -0.077432_wp, &
                      0.154031_wp, 0.153768_wp, 0.161958_wp, 0.163987_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.077233_wp, -0.077241_wp, 0.024499_wp, -0.116775_wp, 0.02447_wp, &
                      0.04994_wp, 0.049948_wp, 0.061195_wp, 0.061197_wp]
         case ("resp")
            charges = [ &
                      -0.19259_wp, -0.193987_wp, -0.079984_wp, -0.136183_wp, -0.084426_wp, &
                      0.16721_wp, 0.16691_wp, 0.175514_wp, 0.177535_wp]
         case ("mbis")
            charges = [ &
                      -0.259636_wp, -0.25978_wp, 0.010882_wp, -0.196172_wp, 0.010726_wp, &
                      0.185847_wp, 0.185878_wp, 0.161114_wp, 0.161139_wp]
         case ("chelpg")
            charges = [ &
                      -0.193007_wp, -0.194517_wp, -0.079815_wp, -0.13621_wp, -0.084347_wp, &
                      0.167497_wp, 0.167224_wp, 0.175563_wp, 0.177611_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.079076_wp, -0.079083_wp, 0.022249_wp, -0.121986_wp, 0.022213_wp, &
                      0.053155_wp, 0.053164_wp, 0.064684_wp, 0.064679_wp]
         case ("resp")
            charges = [ &
                      -0.19935_wp, -0.200751_wp, -0.08393_wp, -0.145763_wp, -0.088431_wp, &
                      0.175116_wp, 0.174813_wp, 0.18313_wp, 0.185165_wp]
         case ("mbis")
            charges = [ &
                      -0.26502_wp, -0.265156_wp, 0.003905_wp, -0.201954_wp, 0.003743_wp, &
                      0.193015_wp, 0.193044_wp, 0.169205_wp, 0.169218_wp]
         case ("chelpg")
            charges = [ &
                      -0.199762_wp, -0.201269_wp, -0.083781_wp, -0.14578_wp, -0.088368_wp, &
                      0.175402_wp, 0.175124_wp, 0.183187_wp, 0.185247_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (177)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.008139_wp, -0.050214_wp, 0.008141_wp, -0.048797_wp, -0.037451_wp, &
                      -0.048796_wp, 0.026863_wp, -0.067893_wp, 0.0312_wp, 0.034124_wp, &
                      0.031201_wp, -0.067894_wp, 0.030606_wp, 0.030613_wp, 0.029467_wp, &
                      0.029467_wp, 0.030607_wp, 0.030615_wp]
         case ("resp")
            charges = [ &
                      0.4515_wp, -0.55618_wp, 0.460654_wp, -0.381633_wp, -0.056531_wp, &
                      -0.379222_wp, 0.197001_wp, -0.49473_wp, 0.169356_wp, 0.134457_wp, &
                      0.168109_wp, -0.465662_wp, 0.121549_wp, 0.1209_wp, 0.123005_wp, &
                      0.130449_wp, 0.128642_wp, 0.128337_wp]
         case ("mbis")
            charges = [ &
                      0.166407_wp, -0.268061_wp, 0.1664_wp, -0.225255_wp, -0.109183_wp, &
                      -0.225252_wp, 0.135589_wp, -0.463075_wp, 0.141988_wp, 0.131445_wp, &
                      0.14199_wp, -0.46308_wp, 0.145656_wp, 0.145674_wp, 0.143711_wp, &
                      0.143712_wp, 0.14566_wp, 0.145672_wp]
         case ("chelpg")
            charges = [ &
                      0.492538_wp, -0.592109_wp, 0.501923_wp, -0.417467_wp, -0.025066_wp, &
                      -0.415014_wp, 0.202633_wp, -0.521899_wp, 0.176442_wp, 0.129491_wp, &
                      0.175205_wp, -0.492486_wp, 0.126775_wp, 0.126154_wp, 0.128884_wp, &
                      0.136375_wp, 0.133968_wp, 0.133652_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.004878_wp, -0.054586_wp, 0.00488_wp, -0.053582_wp, -0.042722_wp, &
                      -0.053581_wp, 0.031115_wp, -0.068875_wp, 0.034783_wp, 0.037576_wp, &
                      0.034784_wp, -0.068876_wp, 0.032595_wp, 0.032602_wp, 0.031904_wp, &
                      0.031905_wp, 0.032595_wp, 0.032604_wp]
         case ("resp")
            charges = [ &
                      0.451357_wp, -0.571994_wp, 0.460531_wp, -0.394139_wp, -0.07231_wp, &
                      -0.391769_wp, 0.210181_wp, -0.511288_wp, 0.180602_wp, 0.1457_wp, &
                      0.179369_wp, -0.48205_wp, 0.128668_wp, 0.12804_wp, 0.130142_wp, &
                      0.137631_wp, 0.135818_wp, 0.135512_wp]
         case ("mbis")
            charges = [ &
                      0.168628_wp, -0.285495_wp, 0.168621_wp, -0.239203_wp, -0.120604_wp, &
                      -0.239199_wp, 0.147399_wp, -0.474211_wp, 0.152087_wp, 0.14127_wp, &
                      0.152089_wp, -0.474216_wp, 0.150822_wp, 0.15084_wp, 0.149752_wp, &
                      0.149753_wp, 0.150826_wp, 0.150838_wp]
         case ("chelpg")
            charges = [ &
                      0.491964_wp, -0.607474_wp, 0.501362_wp, -0.429322_wp, -0.041655_wp, &
                      -0.426913_wp, 0.215738_wp, -0.538514_wp, 0.187568_wp, 0.140923_wp, &
                      0.186346_wp, -0.508944_wp, 0.133923_wp, 0.133323_wp, 0.136059_wp, &
                      0.143593_wp, 0.141169_wp, 0.140853_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.001958_wp, -0.058269_wp, 0.00196_wp, -0.057816_wp, -0.047396_wp, &
                      -0.057814_wp, 0.034918_wp, -0.069566_wp, 0.0378_wp, 0.04039_wp, &
                      0.037802_wp, -0.069567_wp, 0.034324_wp, 0.034331_wp, 0.034142_wp, &
                      0.034143_wp, 0.034324_wp, 0.034333_wp]
         case ("resp")
            charges = [ &
                      0.447746_wp, -0.582536_wp, 0.456949_wp, -0.403329_wp, -0.086932_wp, &
                      -0.400988_wp, 0.221138_wp, -0.523452_wp, 0.189898_wp, 0.155239_wp, &
                      0.188675_wp, -0.494049_wp, 0.134553_wp, 0.133945_wp, 0.136202_wp, &
                      0.143733_wp, 0.141758_wp, 0.141449_wp]
         case ("mbis")
            charges = [ &
                      0.169893_wp, -0.300018_wp, 0.169886_wp, -0.251207_wp, -0.130516_wp, &
                      -0.251202_wp, 0.15768_wp, -0.483005_wp, 0.16067_wp, 0.149483_wp, &
                      0.160672_wp, -0.48301_wp, 0.155136_wp, 0.155155_wp, 0.155043_wp, &
                      0.155045_wp, 0.15514_wp, 0.155154_wp]
         case ("chelpg")
            charges = [ &
                      0.48801_wp, -0.617662_wp, 0.497431_wp, -0.437991_wp, -0.056918_wp, &
                      -0.435612_wp, 0.226641_wp, -0.550723_wp, 0.196773_wp, 0.150606_wp, &
                      0.19556_wp, -0.520997_wp, 0.139828_wp, 0.139248_wp, 0.142148_wp, &
                      0.149722_wp, 0.147128_wp, 0.146808_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (178)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.068819_wp, 0.002915_wp, -0.043971_wp, -0.042897_wp, -0.042897_wp, &
                      -0.043961_wp, 0.002896_wp, -0.068822_wp, 0.029559_wp, 0.029573_wp, &
                      0.030241_wp, 0.030126_wp, 0.033295_wp, 0.033309_wp, 0.030123_wp, &
                      0.029574_wp, 0.030224_wp, 0.029536_wp]
         case ("resp")
            charges = [ &
                      -0.508604_wp, 0.231784_wp, -0.281625_wp, -0.153448_wp, -0.161759_wp, &
                      -0.275148_wp, 0.224067_wp, -0.52494_wp, 0.138645_wp, 0.138229_wp, &
                      0.138298_wp, 0.159494_wp, 0.142188_wp, 0.143877_wp, 0.159654_wp, &
                      0.14254_wp, 0.143473_wp, 0.143275_wp]
         case ("mbis")
            charges = [ &
                      -0.436372_wp, 0.079739_wp, -0.18637_wp, -0.142105_wp, -0.142173_wp, &
                      -0.186296_wp, 0.07975_wp, -0.436279_wp, 0.13922_wp, 0.139149_wp, &
                      0.141301_wp, 0.131981_wp, 0.1335_wp, 0.133534_wp, 0.131963_wp, &
                      0.139091_wp, 0.141232_wp, 0.139137_wp]
         case ("chelpg")
            charges = [ &
                      -0.527079_wp, 0.239129_wp, -0.285026_wp, -0.153444_wp, -0.162392_wp, &
                      -0.27796_wp, 0.230929_wp, -0.543503_wp, 0.142991_wp, 0.142563_wp, &
                      0.143066_wp, 0.160424_wp, 0.142577_wp, 0.144386_wp, 0.160496_wp, &
                      0.146935_wp, 0.148242_wp, 0.147667_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.069785_wp, 0.000338_wp, -0.048968_wp, -0.049047_wp, -0.049058_wp, &
                      -0.048966_wp, 0.000321_wp, -0.069789_wp, 0.032181_wp, 0.032193_wp, &
                      0.032497_wp, 0.034112_wp, 0.036516_wp, 0.036529_wp, 0.034104_wp, &
                      0.032193_wp, 0.032469_wp, 0.032167_wp]
         case ("resp")
            charges = [ &
                      -0.531286_wp, 0.234955_wp, -0.297813_wp, -0.168669_wp, -0.177289_wp, &
                      -0.291178_wp, 0.226723_wp, -0.54735_wp, 0.148183_wp, 0.147756_wp, &
                      0.146342_wp, 0.172268_wp, 0.15345_wp, 0.155216_wp, 0.172437_wp, &
                      0.152007_wp, 0.151468_wp, 0.15278_wp]
         case ("mbis")
            charges = [ &
                      -0.449553_wp, 0.08033_wp, -0.200999_wp, -0.155408_wp, -0.155518_wp, &
                      -0.200929_wp, 0.080348_wp, -0.449461_wp, 0.145855_wp, 0.145779_wp, &
                      0.147308_wp, 0.143229_wp, 0.143525_wp, 0.143564_wp, 0.143206_wp, &
                      0.145721_wp, 0.147219_wp, 0.145789_wp]
         case ("chelpg")
            charges = [ &
                      -0.549926_wp, 0.242338_wp, -0.301183_wp, -0.168715_wp, -0.177935_wp, &
                      -0.293989_wp, 0.233637_wp, -0.566064_wp, 0.152569_wp, 0.152131_wp, &
                      0.151153_wp, 0.17319_wp, 0.153858_wp, 0.155736_wp, 0.173277_wp, &
                      0.156438_wp, 0.156275_wp, 0.157209_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.070313_wp, -0.00157_wp, -0.05276_wp, -0.05395_wp, -0.053971_wp, &
                      -0.052766_wp, -0.001588_wp, -0.070318_wp, 0.034336_wp, 0.034345_wp, &
                      0.034195_wp, 0.037081_wp, 0.038687_wp, 0.038698_wp, 0.037067_wp, &
                      0.034344_wp, 0.034158_wp, 0.034331_wp]
         case ("resp")
            charges = [ &
                      -0.54646_wp, 0.235727_wp, -0.308766_wp, -0.180621_wp, -0.189433_wp, &
                      -0.302056_wp, 0.227096_wp, -0.562282_wp, 0.155406_wp, 0.15497_wp, &
                      0.151993_wp, 0.181519_wp, 0.16159_wp, 0.163404_wp, 0.181704_wp, &
                      0.159161_wp, 0.157073_wp, 0.159976_wp]
         case ("mbis")
            charges = [ &
                      -0.458893_wp, 0.080566_wp, -0.211748_wp, -0.165792_wp, -0.165938_wp, &
                      -0.211686_wp, 0.080591_wp, -0.458803_wp, 0.15102_wp, 0.150936_wp, &
                      0.151669_wp, 0.151586_wp, 0.150743_wp, 0.150787_wp, 0.151556_wp, &
                      0.150876_wp, 0.151564_wp, 0.15097_wp]
         case ("chelpg")
            charges = [ &
                      -0.565208_wp, 0.243133_wp, -0.312113_wp, -0.1807_wp, -0.190089_wp, &
                      -0.304862_wp, 0.23404_wp, -0.581091_wp, 0.159819_wp, 0.159372_wp, &
                      0.156832_wp, 0.182435_wp, 0.16201_wp, 0.163932_wp, 0.182542_wp, &
                      0.163615_wp, 0.161904_wp, 0.164428_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (179)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.003437_wp, -0.043542_wp, -0.043513_wp, 0.003433_wp, -0.043517_wp, &
                      -0.043537_wp, -0.068273_wp, -0.068272_wp, 0.031167_wp, 0.031205_wp, &
                      0.031187_wp, 0.03118_wp, 0.029764_wp, 0.029996_wp, 0.029761_wp, &
                      0.030009_wp, 0.029755_wp, 0.029755_wp]
         case ("resp")
            charges = [ &
                      0.341579_wp, -0.286171_wp, -0.281328_wp, 0.340674_wp, -0.28194_wp, &
                      -0.282076_wp, -0.48763_wp, -0.485822_wp, 0.163596_wp, 0.162093_wp, &
                      0.16213_wp, 0.163109_wp, 0.127596_wp, 0.130848_wp, 0.128229_wp, &
                      0.130515_wp, 0.127023_wp, 0.127575_wp]
         case ("mbis")
            charges = [ &
                      0.134413_wp, -0.190173_wp, -0.19016_wp, 0.134388_wp, -0.190152_wp, &
                      -0.190125_wp, -0.457217_wp, -0.457202_wp, 0.136035_wp, 0.136154_wp, &
                      0.136103_wp, 0.136068_wp, 0.143491_wp, 0.143946_wp, 0.143503_wp, &
                      0.143974_wp, 0.143463_wp, 0.143486_wp]
         case ("chelpg")
            charges = [ &
                      0.35613_wp, -0.291651_wp, -0.286906_wp, 0.354309_wp, -0.286821_wp, &
                      -0.28705_wp, -0.510323_wp, -0.507607_wp, 0.165349_wp, 0.163863_wp, &
                      0.163775_wp, 0.164772_wp, 0.132768_wp, 0.135985_wp, 0.133389_wp, &
                      0.135471_wp, 0.132007_wp, 0.13254_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -2.2e-05_wp, -0.048009_wp, -0.047968_wp, -2.7e-05_wp, -0.047974_wp, &
                      -0.048003_wp, -0.069296_wp, -0.069296_wp, 0.034881_wp, 0.03492_wp, &
                      0.034902_wp, 0.034895_wp, 0.031913_wp, 0.031668_wp, 0.031916_wp, &
                      0.031679_wp, 0.031909_wp, 0.031906_wp]
         case ("resp")
            charges = [ &
                      0.34443_wp, -0.301466_wp, -0.29639_wp, 0.343139_wp, -0.296867_wp, &
                      -0.297124_wp, -0.50531_wp, -0.502744_wp, 0.175486_wp, 0.173893_wp, &
                      0.173893_wp, 0.17495_wp, 0.134717_wp, 0.138026_wp, 0.135365_wp, &
                      0.137514_wp, 0.133967_wp, 0.134522_wp]
         case ("mbis")
            charges = [ &
                      0.13502_wp, -0.203158_wp, -0.203109_wp, 0.134997_wp, -0.203112_wp, &
                      -0.203104_wp, -0.467764_wp, -0.467752_wp, 0.146232_wp, 0.146343_wp, &
                      0.146292_wp, 0.146265_wp, 0.148957_wp, 0.148487_wp, 0.148984_wp, &
                      0.148514_wp, 0.148942_wp, 0.14896_wp]
         case ("chelpg")
            charges = [ &
                      0.359079_wp, -0.306982_wp, -0.302005_wp, 0.356864_wp, -0.30178_wp, &
                      -0.302131_wp, -0.528164_wp, -0.524682_wp, 0.17725_wp, 0.175675_wp, &
                      0.17555_wp, 0.176624_wp, 0.139926_wp, 0.1432_wp, 0.140562_wp, &
                      0.142506_wp, 0.138986_wp, 0.139522_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.003337_wp, -0.052118_wp, -0.052065_wp, -0.003342_wp, -0.052073_wp, &
                      -0.052112_wp, -0.070111_wp, -0.070111_wp, 0.038263_wp, 0.038305_wp, &
                      0.038285_wp, 0.038278_wp, 0.033925_wp, 0.033206_wp, 0.033934_wp, &
                      0.033216_wp, 0.033927_wp, 0.033922_wp]
         case ("resp")
            charges = [ &
                      0.344284_wp, -0.314539_wp, -0.309234_wp, 0.342581_wp, -0.309582_wp, &
                      -0.309942_wp, -0.519267_wp, -0.515923_wp, 0.186101_wp, 0.184425_wp, &
                      0.184391_wp, 0.185516_wp, 0.140969_wp, 0.144326_wp, 0.141627_wp, &
                      0.143633_wp, 0.140037_wp, 0.140597_wp]
         case ("mbis")
            charges = [ &
                      0.134804_wp, -0.214838_wp, -0.214756_wp, 0.134782_wp, -0.214768_wp, &
                      -0.214782_wp, -0.476595_wp, -0.476587_wp, 0.155476_wp, 0.155582_wp, &
                      0.155531_wp, 0.15551_wp, 0.153891_wp, 0.152502_wp, 0.15393_wp, &
                      0.152526_wp, 0.153887_wp, 0.1539_wp]
         case ("chelpg")
            charges = [ &
                      0.359_wp, -0.320079_wp, -0.314875_wp, 0.356366_wp, -0.314515_wp, &
                      -0.31497_wp, -0.542239_wp, -0.537973_wp, 0.187873_wp, 0.186215_wp, &
                      0.186054_wp, 0.187198_wp, 0.146207_wp, 0.149528_wp, 0.146853_wp, &
                      0.148651_wp, 0.145083_wp, 0.145623_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (180)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.201446_wp, 0.084368_wp, -0.064577_wp, -0.031868_wp, -0.053035_wp, &
                      -0.029429_wp, -0.052743_wp, 0.168054_wp, 0.032248_wp, 0.036962_wp, &
                      0.033056_wp, 0.037781_wp, 0.040627_wp]
         case ("resp")
            charges = [ &
                      -0.577481_wp, 0.439659_wp, -0.36921_wp, -0.045569_wp, -0.213837_wp, &
                      -0.0979_wp, -0.261814_wp, 0.403652_wp, 0.158694_wp, 0.125294_wp, &
                      0.136996_wp, 0.131125_wp, 0.170388_wp]
         case ("mbis")
            charges = [ &
                      -0.588398_wp, 0.421647_wp, -0.36244_wp, -0.04844_wp, -0.214496_wp, &
                      -0.070764_wp, -0.297518_wp, 0.441476_wp, 0.150499_wp, 0.127153_wp, &
                      0.142787_wp, 0.130682_wp, 0.167812_wp]
         case ("chelpg")
            charges = [ &
                      -0.582038_wp, 0.462729_wp, -0.392135_wp, -0.025659_wp, -0.232265_wp, &
                      -0.079894_wp, -0.282997_wp, 0.404068_wp, 0.163072_wp, 0.121926_wp, &
                      0.140307_wp, 0.128279_wp, 0.174605_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.210092_wp, 0.080114_wp, -0.064757_wp, -0.035488_wp, -0.057759_wp, &
                      -0.033905_wp, -0.05879_wp, 0.179278_wp, 0.040536_wp, 0.041047_wp, &
                      0.036161_wp, 0.041493_wp, 0.042161_wp]
         case ("resp")
            charges = [ &
                      -0.615293_wp, 0.45795_wp, -0.384006_wp, -0.054053_wp, -0.231485_wp, &
                      -0.103025_wp, -0.286221_wp, 0.434523_wp, 0.176558_wp, 0.135722_wp, &
                      0.147829_wp, 0.141507_wp, 0.179994_wp]
         case ("mbis")
            charges = [ &
                      -0.616511_wp, 0.425096_wp, -0.369626_wp, -0.059217_wp, -0.226812_wp, &
                      -0.078432_wp, -0.314543_wp, 0.46813_wp, 0.167018_wp, 0.137426_wp, &
                      0.151929_wp, 0.140221_wp, 0.175321_wp]
         case ("chelpg")
            charges = [ &
                      -0.619817_wp, 0.480751_wp, -0.406538_wp, -0.034632_wp, -0.249491_wp, &
                      -0.085361_wp, -0.307123_wp, 0.434934_wp, 0.180863_wp, 0.132465_wp, &
                      0.15106_wp, 0.138722_wp, 0.184166_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.217732_wp, 0.076138_wp, -0.063858_wp, -0.03821_wp, -0.062131_wp, &
                      -0.038533_wp, -0.065235_wp, 0.189464_wp, 0.049011_wp, 0.045046_wp, &
                      0.038865_wp, 0.044507_wp, 0.042667_wp]
         case ("resp")
            charges = [ &
                      -0.649073_wp, 0.473298_wp, -0.394044_wp, -0.061162_wp, -0.248113_wp, &
                      -0.108108_wp, -0.31008_wp, 0.46244_wp, 0.193449_wp, 0.145384_wp, &
                      0.157651_wp, 0.150664_wp, 0.187696_wp]
         case ("mbis")
            charges = [ &
                      -0.641044_wp, 0.427551_wp, -0.373572_wp, -0.068491_wp, -0.238142_wp, &
                      -0.085806_wp, -0.331661_wp, 0.491632_wp, 0.182927_wp, 0.146953_wp, &
                      0.16005_wp, 0.148488_wp, 0.181115_wp]
         case ("chelpg")
            charges = [ &
                      -0.65357_wp, 0.495877_wp, -0.416259_wp, -0.042132_wp, -0.265776_wp, &
                      -0.090731_wp, -0.330746_wp, 0.462848_wp, 0.197696_wp, 0.142215_wp, &
                      0.160817_wp, 0.147932_wp, 0.191829_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (181)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (182)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("resp")
            charges = [ &
                      -0.001231_wp, 0.001231_wp]
         case ("mbis")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("chelpg")
            charges = [ &
                      -0.001232_wp, 0.001232_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("resp")
            charges = [ &
                      -0.001274_wp, 0.001274_wp]
         case ("mbis")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("chelpg")
            charges = [ &
                      -0.001275_wp, 0.001275_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("resp")
            charges = [ &
                      -0.001294_wp, 0.001294_wp]
         case ("mbis")
            charges = [ &
                      7e-06_wp, 7e-06_wp]
         case ("chelpg")
            charges = [ &
                      -0.001295_wp, 0.001295_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (183)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.122356_wp, 0.122356_wp]
         case ("resp")
            charges = [ &
                      -0.226794_wp, 0.226794_wp]
         case ("mbis")
            charges = [ &
                      -0.307111_wp, 0.307111_wp]
         case ("chelpg")
            charges = [ &
                      -0.226868_wp, 0.226868_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.145475_wp, 0.145475_wp]
         case ("resp")
            charges = [ &
                      -0.276864_wp, 0.276864_wp]
         case ("mbis")
            charges = [ &
                      -0.353693_wp, 0.353693_wp]
         case ("chelpg")
            charges = [ &
                      -0.276942_wp, 0.276942_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.147541_wp, 0.147541_wp]
         case ("resp")
            charges = [ &
                      -0.281458_wp, 0.281458_wp]
         case ("mbis")
            charges = [ &
                      -0.357792_wp, 0.357792_wp]
         case ("chelpg")
            charges = [ &
                      -0.281537_wp, 0.281537_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (184)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (185)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (186)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      -0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      -0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      -0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      -0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case (187)
      select case (trim(adjustl(to_lower(environment))))
      case ("gas")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("solvent")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case ("conductor")
         select case (trim(adjustl(to_lower(model))))
         case ("hirshfeld")
            charges = [ &
                      0.0_wp]
         case ("resp")
            charges = [ &
                      0.0_wp]
         case ("mbis")
            charges = [ &
                      0.0_wp]
         case ("chelpg")
            charges = [ &
                      0.0_wp]
         case default
            call fatal_error(error, message="Unknown solvent charge model: "//trim(model))
         end select
      case default
         call fatal_error(error, message="Unknown solvent environment: "//trim(environment))
      end select
   case default
      write (msg, "(a,i0)") "No charge data available for solvent ID ", solvent_id
      call fatal_error(error, message=trim(msg))
   end select
end subroutine get_solvent_charges
