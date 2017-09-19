!
! ---------------------------------------------------------------------------------------
!
!                                   Module - Data_DNA
!
!                                                                    Updated : 2017/03/27
!
! Comments: The data structure of DNA data
!
! Script written by Hyungmin Jun (hyungminjun@outlook.com)
! Copyright Hyungmin Jun, 2017. All rights reserved.
!
! ---------------------------------------------------------------------------------------
!
module Data_DNA

! ---------------------------------------------------------------------------------------

    ! BaseType structure
    type :: BaseType
        integer :: id           ! Base ID
        integer :: node         ! Heritage node ID
        integer :: up           ! Upward base ID
        integer :: dn           ! Downward base ID
        integer :: xover        ! Possible crossover ID
        integer :: across       ! Complementary base ID
        integer :: strand       ! Strand ID

        double precision :: pos(3)  ! Position vector
    end type BaseType

! ---------------------------------------------------------------------------------------

    ! TopType structure (base information)
    type :: TopType
        integer   :: id         ! ID
        integer   :: node       ! Node ID
        integer   :: up         ! Upward strand ID
        integer   :: dn         ! Downward strand ID
        integer   :: xover      ! Crossover ID
        integer   :: across     ! Base pair ID
        integer   :: strand     ! Strand ID
        integer   :: address    ! Address number
        logical   :: b_14nt     ! nt of the 14nt seed
        character :: seq        ! Sequence

        double precision :: pos(3)  ! Position vector
    end type TopType

! ---------------------------------------------------------------------------------------

    ! StrandType structure
    type :: StrandType
        integer      :: n_base         ! The number of bases in this strand
        integer      :: n_14nt, n_4nt  ! The number of 14nt and 4nt seeds
        logical      :: b_circular     ! Is it circular strand
        character(4) :: type1          ! Type 1 - Scaffold or staple
        character(6) :: type2          ! Type 2 - Edge or vertex strand

        integer, allocatable :: base(:)
    end type StrandType

! ---------------------------------------------------------------------------------------

    ! DNA data type structure
    type :: DNAType
        integer :: n_base_scaf      ! The number of bases in scaffold strands
        integer :: n_base_stap      ! The number of bases in staple strands
        integer :: n_xover_scaf     ! Possible centered scaffold crossover
        integer :: n_xover_stap     ! Staple crossover
        integer :: n_sxover_stap    ! Single staple crossover
        integer :: n_scaf           ! The number of scaffold strands
        integer :: n_stap           ! The number of staple strands
        integer :: n_top            ! The number of top data
        integer :: n_strand         ! The number of strands

        double precision :: len_ave_stap                ! Average staple length
        integer :: n_14nt, n_s14nt, n_4nt, n_only_4nt   ! The number of 14nt and 4nt seeds
        integer :: n_nt_14nt, n_nt_4nt                  ! The number of nucleotides in 14nt and 4nt domains
        integer :: n_tot_region, n_tot_14nt, n_tot_4nt  ! the total number of 14nt and 4nt seeds
        integer :: len_min_stap, len_max_stap
        integer :: n_unpaired_scaf, n_nt_unpaired_scaf
        integer :: n_unpaired_stap, n_nt_unpaired_stap
        integer :: graph_node, graph_edge
        integer :: min_xover_scaf, min_xover_stap

        type(BaseType),   allocatable :: base_scaf(:)       ! Base in scaffold strand
        type(BaseType),   allocatable :: base_stap(:)       ! Base in staple strand
        type(TopType),    allocatable :: top(:)             ! dnatop
        type(StrandType), allocatable :: strand(:)          ! Strand
        integer,          allocatable :: order_stap(:,:)    ! Staple ordering
    end type DNAType

! ---------------------------------------------------------------------------------------

end module Data_DNA